import requests
from bs4 import BeautifulSoup
import csv
import re
import time

def scrape_gsea(
    browse_url="https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?letter=A",
    output_tsv="gsea_results.tsv",
    sleep_seconds=0.5
):
    """
    A script to:
      1) Fetch a "browse" or "search" page on the GSEA/MSigDB site for human gene sets.
      2) Parse the gene-set table (id="geneSetTable") to get detail-page links.
      3) For each detail page, scrape metadata from the "lists4 human" table.
      4) Download the associated TSV for that gene set, handling both:
         - Row-based TSVs (header includes "GENE_SYMBOL\tGENE_ID"), typically for microRNA sets
         - Key-value TSVs (one line per field, e.g. "GENE_SYMBOLS\tTSPAN1,MPZL2,..."), often for curated sets
      5) Write all results (metadata + gene symbols) to a single TSV file.
    """

    # --------------------------------------------------------
    # 1. Fetch the "browse" (or "search") page that lists gene sets
    # --------------------------------------------------------
    print(f"Fetching browse page: {browse_url}")
    browse_resp = requests.get(browse_url)
    browse_resp.raise_for_status()
    soup = BeautifulSoup(browse_resp.text, "html.parser")

    # The table that contains gene-set links is typically:
    # <table id="geneSetTable" class="lists2 human"> ...
    table = soup.find("table", id="geneSetTable")
    if not table:
        print("Could not find a table with id='geneSetTable'. Aborting.")
        return

    # Collect all <a> links within the table
    detail_links = []
    for td in table.find_all("td"):
        for link_tag in td.find_all("a", href=True):
            relative_url = link_tag["href"]  # e.g. "msigdb/human/geneset/XXX.html"
            # Make a full URL
            full_url = "https://www.gsea-msigdb.org/gsea/" + relative_url
            detail_links.append(full_url)

    print(f"Found {len(detail_links)} gene-set links on the page.")

    # --------------------------------------------------------
    # 2. Prepare result storage and define output columns
    # --------------------------------------------------------
    results = []
    fieldnames = [
        "STANDARD_NAME",
        "SYSTEMATIC_NAME",
        "COLLECTION",
        "MSIGDB_URL",
        "NAMESPACE",
        "DESCRIPTION_BRIEF",
        "DESCRIPTION_FULL",
        "PMID",
        "GEOID",
        "AUTHORS",
        "CONTRIBUTOR",
        "CONTRIBUTOR_ORG",
        "EXACT_SOURCE",
        "FILTERED_BY_SIMILARITY",
        "EXTERNAL_NAMES_FOR_SIMILAR_TERMS",
        "EXTERNAL_DETAILS_URL",
        "SOURCE_MEMBERS",
        "GENE_SYMBOLS",
        "FOUNDER_NAMES"
    ]

    # --------------------------------------------------------
    # 3. Visit each gene-set detail page and gather metadata
    # --------------------------------------------------------
    for i, detail_url in enumerate(detail_links, start=1):
        print(f"[{i}/{len(detail_links)}] Fetching detail page: {detail_url}")
        resp = requests.get(detail_url)
        resp.raise_for_status()
        detail_soup = BeautifulSoup(resp.text, "html.parser")

        # Initialize our metadata dict
        meta = {fn: "" for fn in fieldnames}
        meta["MSIGDB_URL"] = detail_url  # always store the detail page link

        # The detail page typically has a table with class="lists4 human"
        detail_table = detail_soup.find("table", class_="lists4 human")
        if detail_table:
            # Each row has <th> for the label and <td> for the value
            rows = detail_table.find_all("tr")
            for tr in rows:
                th = tr.find("th")
                td = tr.find("td")
                if not th or not td:
                    continue

                label = th.get_text(strip=True).lower()
                value = td.get_text(" ", strip=True)  # join <br> with spaces

                # Simple matching for known labels
                if "standard name" in label:
                    meta["STANDARD_NAME"] = value
                elif "systematic name" in label:
                    meta["SYSTEMATIC_NAME"] = value
                elif "collection" in label:
                    meta["COLLECTION"] = value.replace("\n", " ")
                elif "identifier namespace" in label or "source platform" in label:
                    meta["NAMESPACE"] = value
                elif "brief description" in label:
                    meta["DESCRIPTION_BRIEF"] = value
                elif "full description" in label:
                    meta["DESCRIPTION_FULL"] = value
                elif "source publication" in label:
                    # Could contain Pubmed link or authors
                    pmid_link = td.find("a", href=re.compile("pubmed"))
                    if pmid_link:
                        meta["PMID"] = pmid_link.text.strip().replace("Pubmed","").strip()
                    # Check if "Authors:" is in the text
                    authors_idx = value.lower().find("authors:")
                    if authors_idx != -1:
                        meta["AUTHORS"] = value[authors_idx + len("authors:"):].strip()
                elif "exact source" in label:
                    meta["EXACT_SOURCE"] = value
                elif "filtered by similarity" in label:
                    meta["FILTERED_BY_SIMILARITY"] = value
                elif "external links" in label:
                    meta["EXTERNAL_DETAILS_URL"] = value
                elif "contributed by" in label:
                    # Typically "Name (Org)"
                    # e.g. "Dharmesh D. Bhuva (Walter and Eliza Hall Institute of Medical Research)"
                    match_paren = re.search(r"(.*)\((.*)\)", value)
                    if match_paren:
                        meta["CONTRIBUTOR"] = match_paren.group(1).strip()
                        meta["CONTRIBUTOR_ORG"] = match_paren.group(2).strip()
                    else:
                        # fallback if no parentheses
                        meta["CONTRIBUTOR"] = value
                elif "founder" in label:
                    meta["FOUNDER_NAMES"] = value
                # etc. Additional fields can be mapped as needed.
        else:
            print("  WARNING: No table with class='lists4 human' found on detail page. Skipping metadata extraction.")

        # --------------------------------------------------------
        # 4. Download the TSV that contains gene info
        # --------------------------------------------------------
        tsv_link_tag = detail_soup.find("a", href=re.compile("download_geneset.jsp.*fileType=TSV"))
        if not tsv_link_tag:
            print("  WARNING: No TSV link found. Cannot collect gene members.")
            results.append(meta)
            time.sleep(sleep_seconds)
            continue

        tsv_url = "https://www.gsea-msigdb.org/gsea/" + tsv_link_tag["href"]
        print(f"  Downloading TSV: {tsv_url}")
        tsv_resp = requests.get(tsv_url)
        tsv_resp.raise_for_status()

        # We'll parse the TSV in two ways:
        #  - If it's row-based with "GENE_SYMBOL\tGENE_ID" as the first line
        #  - Or if it's key-value lines (like "GENE_SYMBOLS\tTSPAN1,MPZL2,...")

        lines = tsv_resp.text.strip().split("\n")
        if not lines:
            print("  WARNING: TSV is empty.")
            results.append(meta)
            time.sleep(sleep_seconds)
            continue

        # Check the first line to see which format we have
        header_line = lines[0].lower()

        # Case A: row-based format
        # e.g. first line might be "GENE_SYMBOL   GENE_ID   ..."
        if "gene_symbol" in header_line and "gene_id" in header_line:
            # We can parse each subsequent line as one gene per row
            gene_symbols = []
            gene_ids = []
            for idx, line in enumerate(lines):
                # skip the header
                if idx == 0:
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    symbol = parts[0].strip()
                    gid = parts[1].strip()
                    gene_symbols.append(symbol)
                    gene_ids.append(gid)
            meta["GENE_SYMBOLS"] = ",".join(gene_symbols)
            meta["SOURCE_MEMBERS"] = ",".join(gene_ids)

        else:
            # Case B: key-value format
            # e.g. lines like:
            #   GENE_SYMBOLS    TSPAN1,MPZL2,VAV3,...
            #   SOURCE_MEMBERS  10103,10205,...
            for line in lines:
                if not line.strip():
                    continue
                kv_parts = line.split("\t", 1)
                if len(kv_parts) != 2:
                    continue
                key = kv_parts[0].strip()
                val = kv_parts[1].strip()

                # We specifically want the lines "GENE_SYMBOLS" and "SOURCE_MEMBERS"
                if key == "GENE_SYMBOLS":
                    meta["GENE_SYMBOLS"] = val
                elif key == "SOURCE_MEMBERS":
                    meta["SOURCE_MEMBERS"] = val
                # If there are other fields in this key-value TSV you want, parse them too.

        # Save the meta for this gene set
        results.append(meta)
        time.sleep(sleep_seconds)

    # --------------------------------------------------------
    # 5. Write everything to a TSV file
    # --------------------------------------------------------
    with open(output_tsv, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for entry in results:
            writer.writerow(entry)

    print(f"\nDone! Collected {len(results)} gene sets. Results in: {output_tsv}")


if __name__ == "__main__":
    """
    Example usage:

    1) Browse gene sets starting with letter A:
       python gsea.py

    2) Or pass a different URL, e.g. to search for 'TGF':
       scrape_gsea(
         browse_url="https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?geneSetName=TGF&Search=Search",
         output_tsv="gsea_tgf_results.tsv"
       )

    Adjust 'sleep_seconds' to be polite. 
    """
    # Default: scrape sets starting with letter A
    scrape_gsea()