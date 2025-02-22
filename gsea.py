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
    Example: This function starts from a "browse" page listing gene sets (like 'letter=A'),
    parses the table with id="geneSetTable", and follows each gene set link to its detail page.
    On the detail page, we extract key metadata and the link for the 'TSV' download.
    Finally, we write the results into a .tsv file with columns that match your example.
    
    You can adapt 'browse_url' to something like:
      https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?geneSetName=TGF&Search=Search
    if you prefer searching by a specific geneSetName or partial name.
    """
    # -------------------------------------------------------------------
    # 1. Fetch the browse page with the geneSet table.
    # -------------------------------------------------------------------
    print(f"Fetching browse page: {browse_url}")
    browse_resp = requests.get(browse_url)
    browse_resp.raise_for_status()
    soup = BeautifulSoup(browse_resp.text, "html.parser")

    # -------------------------------------------------------------------
    # 2. Locate the table with id="geneSetTable".
    #    Each <td> can contain multiple <a> tags, each referencing a gene set detail page.
    # -------------------------------------------------------------------
    table = soup.find("table", id="geneSetTable")
    if not table:
        print("Could not find a table with id='geneSetTable'. Aborting.")
        return

    # Collect all gene-set detail links from the table
    detail_links = []
    all_tds = table.find_all("td")
    for td in all_tds:
        for link_tag in td.find_all("a", href=True):
            # e.g. href="msigdb/human/geneset/XXXX.html"
            relative_url = link_tag["href"]
            # Build the absolute URL
            full_url = "https://www.gsea-msigdb.org/gsea/" + relative_url
            detail_links.append(full_url)

    print(f"Found {len(detail_links)} gene-set links on the browse page.")

    # -------------------------------------------------------------------
    # 3. Prepare our list of results, and the output columns
    # -------------------------------------------------------------------
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

    # -------------------------------------------------------------------
    # 4. Visit each gene-set detail link and scrape data
    # -------------------------------------------------------------------
    for i, detail_url in enumerate(detail_links, start=1):
        print(f"[{i}/{len(detail_links)}] Fetching detail page: {detail_url}")
        resp = requests.get(detail_url)
        resp.raise_for_status()
        detail_soup = BeautifulSoup(resp.text, "html.parser")

        # We'll store all our metadata in a dictionary
        meta = {fn: "" for fn in fieldnames}

        # 4a. Basic fields
        meta["MSIGDB_URL"] = detail_url

        # The detail page has a table with class="lists4 human" containing <th> + <td>
        detail_table = detail_soup.find("table", class_="lists4 human")
        if not detail_table:
            # Some sets might not follow this pattern, or the site changed:
            print(" WARNING: No 'lists4 human' table found. Skipping metadata extraction.")
            # We'll skip or store empty fields, but continue to attempt TSV link anyway.
        else:
            # Each row <tr> typically has <th> for the label and <td> for the value
            rows = detail_table.find_all("tr")
            for tr in rows:
                th = tr.find("th")
                td = tr.find("td")
                if not th or not td:
                    continue
                label = th.get_text(strip=True).lower()
                value = td.get_text(" ", strip=True)  # joined by spaces

                # Match each label to a field in meta
                if "standard name" in label:
                    meta["STANDARD_NAME"] = value
                elif "systematic name" in label:
                    meta["SYSTEMATIC_NAME"] = value
                elif "collection" in label:
                    # e.g.: "C2: Curated  CGP: Chemical and Genetic Perturbations"
                    meta["COLLECTION"] = value.replace("\n", " ")
                elif "source platform" in label or "identifier namespace" in label:
                    meta["NAMESPACE"] = value
                elif "brief description" in label:
                    meta["DESCRIPTION_BRIEF"] = value
                elif "full description" in label:
                    meta["DESCRIPTION_FULL"] = value
                elif "source publication" in label:
                    # May contain "Pubmed 1234" or "Authors: X,Y,Z"
                    # We might parse out PMID etc. For a quick approach, just store the raw text.
                    # If there's a link e.g. <a href="...">Pubmed 28119430</a> we can parse it
                    pmid_link = td.find("a", href=re.compile("pubmed"))
                    if pmid_link:
                        meta["PMID"] = pmid_link.text.strip().replace("Pubmed","").strip()
                    # authors sometimes in the same row
                    # see if there's a string like "Authors: Foroutan M,Cursons J, etc."
                    # We'll just store them in meta["AUTHORS"] if found
                    authors_idx = value.lower().find("authors:")
                    if authors_idx != -1:
                        meta["AUTHORS"] = value[authors_idx + len("authors:"):].strip()
                elif "exact source" in label:
                    meta["EXACT_SOURCE"] = value
                elif "filtered by similarity" in label:
                    meta["FILTERED_BY_SIMILARITY"] = value
                elif "external links" in label:
                    # This might contain clickable links
                    meta["EXTERNAL_DETAILS_URL"] = value
                elif "contributed by" in label:
                    # Usually "Name (Institution)" or something similar
                    # e.g. "Dharmesh Bhuva (Walter and Eliza Hall Institute...)"
                    meta["CONTRIBUTOR"] = value
                    # We can try splitting at parentheses if we want to separate contributor from org
                    # But for a quick approach:
                    # "Name (Org)" => 
                    #  contributor = "Name"
                    #  contributor org = "Org"
                    match_paren = re.search(r"(.*)\s+\((.*)\)", value)
                    if match_paren:
                        meta["CONTRIBUTOR"] = match_paren.group(1).strip()
                        meta["CONTRIBUTOR_ORG"] = match_paren.group(2).strip()
                elif "related gene sets" in label:
                    # Not typically in your final TSV, so we skip or parse if you want
                    pass
                elif "founder" in label:
                    meta["FOUNDER_NAMES"] = value
                # etc. Add more elif statements as you see more fields.

        # 4b. Find the link to download the “TSV” for gene members.
        #     The detail page row with <th> “Download gene set” might have <a> for TSV, etc.
        tsv_a_tag = detail_soup.find("a", href=re.compile("download_geneset.jsp.*fileType=TSV"))
        if tsv_a_tag:
            tsv_url = "https://www.gsea-msigdb.org/gsea/" + tsv_a_tag["href"]
            print(f"  Downloading TSV metadata from {tsv_url}")
            tsv_resp = requests.get(tsv_url)
            tsv_resp.raise_for_status()

            gene_ids = []
            gene_symbols = []
            # The TGF example suggests the first line is a header row
            # The columns might be GENE_SYMBOL, GENE_ID, ...
            for line in tsv_resp.text.splitlines():
                if line.lower().startswith("gene_symbol") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    symbol = parts[0].strip()
                    gid = parts[1].strip()
                    gene_symbols.append(symbol)
                    gene_ids.append(gid)
                elif len(parts) == 1:
                    # If there's only a single column with a symbol
                    gene_symbols.append(parts[0].strip())
                # else: handle other columns if needed

            meta["GENE_SYMBOLS"] = ",".join(gene_symbols)
            meta["SOURCE_MEMBERS"] = ",".join(gene_ids)
        else:
            print("  WARNING: No TSV link found on detail page.")

        results.append(meta)
        time.sleep(sleep_seconds)

    # -------------------------------------------------------------------
    # 5. Write all results to a TSV file in the desired column order
    # -------------------------------------------------------------------
    with open(output_tsv, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for entry in results:
            writer.writerow(entry)

    print(f"Done! Collected {len(results)} gene sets. Output saved to: {output_tsv}")


if __name__ == "__main__":
    # Example usage:
    # 1) By letter:
    #    python script.py => scrapes the "letter=A" browse
    # 2) By search term:
    #    scrape_gsea("https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?geneSetName=TGF&Search=Search")
    #
    scrape_gsea()