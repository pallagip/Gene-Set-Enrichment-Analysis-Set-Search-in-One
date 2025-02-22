# Gene-Set-Enrichment-Analysis-Set-Search-in-One
A Python web scraper that automates collecting gene set metadata (and associated genes) from the GSEA/MSigDB website into a consolidated TSV file.

This script is a Python-based web scraper designed to streamline the collection of gene set information from the GSEA/MSigDB website. By leveraging the BeautifulSoup and requests libraries, it identifies all relevant gene set links on a “browse” or “search” page, then visits each detail page to extract fields such as standard name, systematic name, and descriptive metadata. Additionally, it fetches the TSV file containing source members and gene symbols, consolidating everything into a convenient TSV output file.

The intention behind this code is to automate a process that would otherwise require manual navigation and copying of data from dozens or even hundreds of GSEA pages. Researchers focusing on certain pathways, like TGFβ-related gene sets, can use this script to quickly assemble all metadata, descriptions, and gene listings for further analyses. By reducing manual effort, it not only saves time but also minimizes the chance of transcription errors.
