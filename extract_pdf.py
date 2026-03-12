#!/usr/bin/env python3
"""Extract text from the MATLAB Finite Element Analysis PDF"""

import pypdf

def extract_text_from_pdf(pdf_path, output_path):
    """Extract all text from PDF and save to output file"""
    reader = pypdf.PdfReader(pdf_path)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("MATLAB Codes for Finite Element Analysis - Extracted Text\n")
        f.write("=" * 80 + "\n\n")
        
        for page_num, page in enumerate(reader.pages, 1):
            f.write(f"\n{'='*80}\n")
            f.write(f"PAGE {page_num}\n")
            f.write(f"{'='*80}\n\n")
            
            text = page.extract_text()
            f.write(text)
            f.write("\n\n")
            
    print(f"Text extracted successfully to {output_path}")
    print(f"Total pages: {len(reader.pages)}")

if __name__ == "__main__":
    pdf_path = "2009 - MATLAB Codes for Finite Element Analysis.pdf"
    output_path = "pdf_content.txt"
    extract_text_from_pdf(pdf_path, output_path)
