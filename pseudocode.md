Problem: 

In many cases, amplifying RNA with PCR is necessary to create enough molecules for sequencing. This is not an issue when the goal of the experiment is to view only the sequence of the RNA and not relative abundance; however, in cases where biological number of reads matters, PCR presents a challenge because reads may have the same sequence for a variety of reasons. Reads could have the same sequence because they're duplicates from being replicated, or because they represent a highly expressed RNA within a cell. 

A marker called a Unique Molecular Index (UMI) can be helpful in determining which strand each PCR read comes from. A different UMI can be applied to each molecule and when deduplication is done, reads with the same chromosome, position, strand information, and UMI will be considered duplicates.   

def process_sam_file(file_path: str) -> None:
    known_umis = load_known_umis()  # Assume we have a list of known UMIs
    seen_reads = set()
    
    with open(file_path, 'r') as sam_file, open('output.sam', 'w') as output_file:
        for line in sam_file:
            if line.startswith('@'):  # Header lines
                write_read_to_output(line, output_file)
                continue

            qname, chrom, pos, cigar, seq, qual = parse_sam_line(line)
            umi = extract_umi(qname)

            if not is_valid_umi(umi, known_umis):
                continue  # Skip invalid UMI reads

            read_hash = generate_read_hash(chrom, pos, umi, cigar)
            
            if read_hash not in seen_reads:
                write_read_to_output(line, output_file)
                seen_reads.add(read_hash)

def extract_umi(qname: str) -> str:
    return qname.split(':')[-1]

def is_valid_umi(umi: str, known_umis: list) -> bool:
    return umi in known_umis

def generate_read_hash(chrom: str, pos: int, umi: str, cigar: str) -> str:
    return f'{chrom}_{pos}_{umi}_{cigar}'

def write_read_to_output(read: str, output_file) -> None:
    output_file.write(read)
