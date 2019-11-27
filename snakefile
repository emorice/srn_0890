import json
import pandas as pd
import Bio.SeqIO as SeqIO

rule all:
    input:
        'matches.cmc.results.txt'

rule rep_seqs:
    output:
        'prok_representative_genomes.txt'
    shell:
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_representative_genomes.txt
        """

checkpoint choose_genomes:
    input:
        'prok_representative_genomes.txt'
    output:
        'seqlist.json'
    run:
        reps = pd.read_table(input[0], sep='\t')
        chrs = list(
            reps[reps['#Species/genus'].str.startswith('Staphylococcus')]
            ['Chromosome RefSeq']
            .dropna()
            )
        with open(output[0], 'w') as fd:
            json.dump(chrs, fd)

def input_accessions(wildcards):
    with open(checkpoints.choose_genomes.get().output[0]) as fd:
        acc_list = json.load(fd)
    return [
        '%s.psl' % acc for acc in acc_list
        ]

rule fetch_seqs:
    output:
        "{acc}.fa"
    shell:
        """
        wget 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={wildcards.acc}&rettype=fasta&retmode=text' -O {output}
        """

rule collect_blat:
    input:
        input_accessions,
        'prok_representative_genomes.txt'
    output:
        'matches.json',
        'matches.fasta'
    run:
        matches = []
        with open(output[-1], 'w') as out_fasta:
            for psl_path in input[:-1]:
                try:
                    seq_matches = (pd
                        .read_table(psl_path, sep='\t', header=None)
                        .rename({
                            8: 'strand',
                            10: 'qsize',
                            11: 'qbegin',
                            12: 'qend',
                            13: 'accession',
                            15: 'begin',
                            16: 'end'
                            }, axis=1)
                        )
                    matches.append(seq_matches)
                    genome = SeqIO.read('%s.fa' % seq_matches['accession'].iloc[0], 'fasta')
                    for _, m in seq_matches.iterrows():
                        seq = genome.seq[
                                (m['begin'] - m['qbegin'] + 1 - 20)
                                :
                                (m['end'] + m['qsize'] - m['qend'])
                                ]
                        if m.strand == '-':
                            seq = seq.reverse_complement()
                        print('>', m['accession'], file=out_fasta)
                        print(seq, file=out_fasta)
                except pd.errors.EmptyDataError:
                    pass
            reps = pd.read_table(input[-1], sep='\t').rename(
                {'Chromosome RefSeq': 'accession'},
                axis=1)
            (pd
                .concat(matches)
                .merge(reps, on='accession')
                .to_csv(output[0])
                )
            
rule fetch_blat:
    output:
        'blat'
    shell:
        """
        wget 'https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat'
        chmod +x blat
        """
rule blat_rna:
    input:
        'blat',
        '{acc}.fa'
    output:
        '{acc}.psl'
    shell:
        """
        ./blat {input[1]} srn_0890.fa {output} -t=dna -q=rna -noHead
        """

rule fetch_infernal:
    output:
        'infernal-1.1.3-linux-intel-gcc/README'
    shell:
        """
        wget 'http://eddylab.org/infernal/infernal-1.1.3-linux-intel-gcc.tar.gz' -O - | tar xvzf -
        """

rule cmpress:
    input:
        'infernal-1.1.3-linux-intel-gcc/README',
        '{cm}',
    output:
        '{cm}.i1f',
        '{cm}.i1m',
        '{cm}.i1i',
        '{cm}.i1p'
    shell:
        """
        infernal-1.1.3-linux-intel-gcc/binaries/cmpress {wildcards.cm}
        """

rule cmcalibrate:
    input:
        '{cm}',
        'infernal-1.1.3-linux-intel-gcc/README',
    output:
        '{cm}c'
    shell:
        """
        cp {input[0]} {input[0]}.bak
        infernal-1.1.3-linux-intel-gcc/binaries/cmcalibrate {input[0]}
        mv {input[0]} {output}
        mv {input[0]}.bak {input[0]}
        """

rule cmsearch:
    input:
        'blast.fa.nhr',
        'infernal-1.1.3-linux-intel-gcc/README',
        '{cm}.i1f'
    output:
        '{cm}.results.txt',
        '{cm}.results.msa',
        '{cm}.results.matches'
    shell:
        """
        infernal-1.1.3-linux-intel-gcc/binaries/cmsearch \
                -o {output[0]} -A {output[1]} --tblout {output[2]} {wildcards.cm} blast.fa 
        """

rule fetch_viennarna:
    output:
        'ViennaRNA-2.4.14/README.md'
    shell:
        """
        wget 'https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz' -O - | tar xzf -
        mkdir vienna
        cd ViennaRNA-2.4.14
        ./configure --prefix=$(pwd)/../vienna
        make
        make install
        """

rule fetch_locarna:
    input:
        'ViennaRNA-2.4.14/README.md'
    output:
        'locarna-1.9.2.1/README'
    shell:
        """
        wget 'https://github.com/s-will/LocARNA/releases/download/v1.9.2.1/locarna-1.9.2.1.tar.gz' -O - | tar xzvf -
	mkdir locarna
        cd locarna-1.9.2.1
        ./configure --with-vrna=../vienna --prefix=$(pwd)/../locarna
        make
	make install
        """

rule first_align:
    input:
        'locarna/bin/mlocarna',
        'matches.fasta'
    output:
        'matches.locarna/results/result.stk'
    shell:
        """
        {input[0]} --stockholm {input[1]} --tgtdir matches.locarna
        """
rule cmbuild:
    input:
        'infernal-1.1.3-linux-intel-gcc/README',
        'matches.locarna/results/result.stk'
    output:
        'matches.cm'
    shell:
        """
        infernal-1.1.3-linux-intel-gcc/binaries/cmbuild {output} 'matches.locarna/results/result.stk'
        """
