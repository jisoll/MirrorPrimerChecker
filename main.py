import os
import argparse
import pandas as pd
import re
from datetime import datetime

def generate_primer_combinations(primer):
    base_map = {
        'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
        'R': 'AG', 'Y': 'CT', 'S': 'GC',
        'W': 'AT', 'K': 'GT', 'M': 'AC',
        'B': 'CGT', 'D': 'AGT', 'H': 'ACT',
        'V': 'ACG', 'N': 'ACGT'
    }
    
    def expand(sequence):
        if sequence:
            first, *rest = sequence
            for char in base_map.get(first, first):
                for item in expand(rest):
                    yield char + item
        else:
            yield ""

    return list(expand(primer.upper().replace(" ", "")))

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                  'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
                  'H': 'D', 'V': 'B', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def process_database(db_path, forward_primers, reverse_primers=None):
    df = pd.read_csv(db_path, sep='\t', header=0)
    results = []

    for _, row in df.iterrows():
        sequence = row['Sequence'].upper()
        matched_forward_primers = [fp for fp in forward_primers if re.search(fp, sequence)]
        matched_reverse_primers = [rp for rp in reverse_primers if re.search(rp, sequence)] if reverse_primers else []
        
        if matched_forward_primers and (not reverse_primers or matched_reverse_primers):
            fwd_match = re.search(matched_forward_primers[0], sequence)
            rev_match = re.search(matched_reverse_primers[0], sequence) if matched_reverse_primers else None
            amplicon_sequence = sequence[fwd_match.start():rev_match.end()] if fwd_match and rev_match else 'N/A'
            amplicon_length = len(amplicon_sequence) if amplicon_sequence != 'N/A' else 0
                         
            results.append({
                'Accession': row['Accession'],
                'GenBankAccession': row['GenBankAccession'],
                'OperonNumber': row['OperonNumber'],
                'Domain': row['Domain'],
                'Phylum': row['Phylum'],
                'Class': row['Class'],
                'Order': row['Order'],
                'Family': row['Family'],
                'Genus': row['Genus'],
                'Species': row['Species'],
                'Amplicon_Sequence': amplicon_sequence,
                'Amplicon_Length': amplicon_length,
                'Forward_Primer': ','.join(matched_forward_primers),
                'Reverse_Primer': ','.join(matched_reverse_primers)
            })
    return results, df

def ensure_directory(directory):
    os.makedirs(directory, exist_ok=True)

def save_results(forward_primers, reverse_primers, results, df, output_folder):
    # Save query primer sequences
    query_primer_file = os.path.join(output_folder, '2_Query_primer_sequences.txt')
    with open(query_primer_file, 'w') as file:
        file.write('# Forward Primer Sequence(s)\n')
        for primer in forward_primers:
            file.write(primer + '\n')
        if reverse_primers:
            file.write('\n# Reverse Primer Sequence(s)\n')
            for primer in reverse_primers:
                file.write(primer + '\n')
    print(f"{query_primer_file} has been saved.")
    
    # Save matched sequences
    sequence_file = os.path.join(output_folder, '3_Database_matches.tsv')
    pd.DataFrame(results).to_csv(sequence_file, index=False, sep='\t')
    print(f"{sequence_file} has been saved.")
    
    # Save summary
    summary_file = os.path.join(output_folder, '1_Summary.txt')
    with open(summary_file, 'w') as f:
        f.write("# Input Primers:\n")
        f.write(f"\t- Forward Primer: 5'-> 3' {forward_primers[0]}\n")
        if reverse_primers:
            f.write(f"\t- Reverse Primer: 5'-> 3' {reverse_primers[0]}\n\n")
        
        if results:
            total_sequences = len(results)
            unique_genomes = pd.DataFrame(results)['Accession'].nunique()
            unique_species = pd.DataFrame(results)['Species'].nunique()
            
            domain_counts = pd.DataFrame(results)['Domain'].value_counts()
            domain_total = domain_counts.sum()
            f.write("# Domain Distribution:\n")
            for domain, count in domain_counts.items():
                percentage = (count / domain_total) * 100
                f.write(f"\t- {domain}: {count} ({percentage:.2f}%)\n")
            f.write("\n")
            
            f.write("# Operon Sequences:\n")
            f.write(f"\t- Percentage of Matches: {total_sequences / len(df) * 100:.2f}%\n")
            f.write(f"\t- Matched: {total_sequences}\n")
            f.write(f"\t- Total Covered: {len(df)}\n\n")
            f.write("# Genomes:\n")
            f.write(f"\t- Matched: {unique_genomes}\n")
            f.write(f"\t- Total Covered: {df['Accession'].nunique()}\n\n")
            f.write("# Species:\n")
            f.write(f"\t- Matched: {unique_species}\n")
            f.write(f"\t- Total Covered: {df['Species'].nunique()}\n\n")
        
    print(f"{summary_file} has been saved.")

def main():
    parser = argparse.ArgumentParser(description="Primer Matching Tool")
    parser.add_argument("forward_primer", help="Forward primer sequence")
    parser.add_argument("reverse_primer", help="Reverse primer sequence")
    parser.add_argument("--db_path", required=True, help="Path to the primer database (TSV format)")
    parser.add_argument("--output_dir", default="output", help="Directory to save results")
    args = parser.parse_args()

    start_time = datetime.now()
    print(f"MirrorPrimerChecker started at: {start_time}\n")

    forward_primers = generate_primer_combinations(args.forward_primer)
    reverse_primers = generate_primer_combinations(reverse_complement(args.reverse_primer)) if args.reverse_primer else None
    
    timestamp = start_time.strftime("%Y%m%d_%H%M%S")
    output_folder = os.path.join(args.output_dir, f"output_{timestamp}")
    ensure_directory(output_folder)

    results, df = process_database(args.db_path, forward_primers, reverse_primers)
    save_results(forward_primers, reverse_primers, results, df, output_folder)

    end_time = datetime.now()
    print(f"\n")
    print(f"MirrorPrimerChecker ended at: {end_time}")
    print(f"Total execution time: {end_time - start_time}\n")

if __name__ == "__main__":
    main()
