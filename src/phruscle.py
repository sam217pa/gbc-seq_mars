# -*- coding: utf-8 -*-

# 1. lancer phred sur la séquence abi en entrée.
# 2. lancer muscle avec pour entrée le fichier d'alignement fst entre
#    le gène synthétique et la séquence sauvage d'une part, et la
#    séquence dans le fichier seq d'autre part.
# 3. analyser la sortie de phred .phd pour extraire l'information de
#    la qualité de la base.

import os.path
import io
import subprocess as sub  # permet de lancer phred et muscle
import click  # construit la cli
from Bio import AlignIO  # permet d'analyser les sorties de muscle.
from Bio import SeqIO  # permet de tester les entrées
from Bio.Seq import Seq
import numpy as np
import pandas as pd

# TODO: définir les chemins de fichiers avec constance.


@click.group()
def phruscle():
    """This program gather a set of function to run phred on a abi file, to run muscle on the
    called bases, and to give a nicely formatted output as csv.

    """


def is_abi(sequence):
    """Check if sequence if a spectrogram from abi file.

    :param sequence: sequence file name. str.

    """
    if os.path.splitext(sequence)[1].lower().endswith(('abi', 'ab1')):
        return True
    else:
        return False


def is_fasta(sequence):
    """Check if sequence is a fasta file.

    :param sequence: name of the fasta file. str

    """
    ext = os.path.splitext(sequence)[1].lower()
    if ext.endswith(('.seq', '.fasta', '.aln', '.fst', '.rev')):
        return True
    else:
        return False


##-----------------------------------------------------------------------------
##                                          BASECALL
##-----------------------------------------------------------------------------
@phruscle.command('basecall', short_help='Make basecall and alignment.')
@click.option('-i',
              '--input',
              type=click.File('rb'),
              default='-',
              help="The abi sequence to base call with phred")
@click.option('-r',
              '--ref',
              type=click.File('r'),
              default='-',
              help="The reference alignment in fasta format")
@click.option('-c',
              '--cutoff',
              default=0.05,
              help="Probability of error of basecalling")
@click.option('--reverse',
              is_flag=True,
              help="Reverse complement sequencing file.")
@click.option('--clustalw',
              is_flag=True,
              help="Muscle output in human readable clustalw format")
def run_phruscle(input, ref, cutoff, reverse, clustalw):
    """This programme uses phred to call bases on the `--input` file to create a fasta file.
    Muscle align this fasta file to the `--ref`, without disrupting the reference,
    to get the localisation of SNPs and their polarity.

    :param input: the abi file to run phred on. str
    :param ref: the reference alignment. str
    :param cutoff: the maximum phred probability of error.
    :param reverse: bool. compute reverse complement of input.
    :param clustalw: bool. muscle output as human readable format clustalw.
    """

    def run_phred(abi_file, cutoff):
        """Run phred on the abi_file, with a probability of base calling error of cutoff.

        :param abi_file: the sanger sequencing file to base call. str
        :param cutoff: probability of error. float.
        """
        # print abi_file.name
        assert is_abi(abi_file.name), "Abi_File is not, well…, an abi file."

        # phred reads and process file listed in the file abi_file via `-if`.
        # we must then create a temporary file containing only the name of the abi_file.
        with open('tempfile', 'w') as tempfile:
            tempfile.write(abi_file.name + "\n")

        # subprocess.call takes a list of argument as abi_file.
        # the following list describes the phred list of argument.
        phred = [
            "phred",
            "-st",
            "fasta",  # sequence type output (default = fasta)
            "-trim_alt",
            "\"\"",  # quality trim
            "-trim_cutoff",
            str(cutoff),  # error probability of trimming
            "-trim_fasta",  # trim sortie fasta.
            "-trim_phd",  # trim sortie phd
            "-s",  # write seq file, append ".seq" to the names of the abi_file files
            "-d",  # write a poly file
            "-p",  # write a phd file
            "-if",
            "tempfile"  # abi_file file
        ]
        phred_call = sub.call(phred)
        os.remove("tempfile")  # supprime le fichier temporaire
        return phred_call

    def run_muscle(phred_output, reference, clw=False):
        """Run muscle on the phred_output sequence, and align it to the reference sequence.
        :param phred_output: the fasta sequence to align.
        :param reference: the reference fasta sequence.

        """
        assert is_fasta(
            phred_output), "La séquence n'est pas un fichier fasta."
        assert is_fasta(
            reference.name), "La référence n'est pas un fichier fasta."

        muscle = [
            "muscle",
            "-quiet",  # non verbose
            "-profile",  # do not disrupt the profile alignment
            # "-clw",
            # "-objscore",
            # "ps",  # matrice de scoring
            # "-maxmb",
            # str(100),  # 100M as max memory
            "-in1",
            phred_output,
            "-in2",
            reference.name,  # profile is reference
            "-out",
            phred_output + ".aln"  # append .aln to phred_output name
        ]
        if clw:
            muscle.append("-clw")

        muscle_call = sub.call(muscle)  # call muscle list on input
        return muscle_call  # return 0 or 1 if muscle is successful.

    def reverse_complement(seq_in, seq_out):
        """Read seq_in and writes its reverse complement to seq_out.

        :param seq_in: sequence to reverse comp. str
        :param seq_out: sequence to write revcomp to. str
        :returns: None.
        :rtype: NoneType.

        """
        record = SeqIO.read(open(seq_in), "fasta")
        SeqIO.write(record.reverse_complement(), open(seq_out, 'w'), 'fasta')

    phred_code = run_phred(input, cutoff)

    phred_output = os.path.basename(input.name) + ".seq"
    if reverse:
        reverse_complement(phred_output, phred_output + ".rev")
        phred_output += ".rev"

    muscle_code = run_muscle(phred_output, ref, clustalw)

    if phred_code == 0 and muscle_code == 0:
        print "Phruscle ran fine."
    elif muscle_code != 0:
        print "Muscle failed."
    else:
        print "Phred failed."

##-----------------------------------------------------------------------------
## TABLE
##-----------------------------------------------------------------------------
# for line in reversed(open("filename").readlines()):
#     print line.rstrip()


@phruscle.command('table',
                  short_help="Make a table of SNP quality and polarity")
@click.option('-i', '--input', help="Input alignment in fasta format.")
@click.option(
    '-p',
    '--phd',
    help="Input sequence quality and position file, output by phred.")
@click.option('-o',
              '--output',
              type=click.File('wb'),
              default='-', # default is print to stdin
              help="Output file in csv format. Default is STDOUT.")
@click.option(
    '-r',
    '--reverse',
    is_flag=True,
    help=
    "Reverse the phd file. Useful when reference and sequence are not in the same direction."
)

def clean_output(input, phd, output, reverse=False):
    """This programme takes a fasta alignment between the experimental sequence and the
    reference alignment, find the quality of the base call in the phd file, and outputs it
    as csv.

    \b
    To read the consensus (cons) output:
    - `.` : the three bases are the same.
    - `x` : the sequenced base matches to the snp base.
    - `X` : the sequenced base matches the WT base.
    - `-` : the sequenced base is misaligned.
    - `N` : the base was not called.

    """

    def parse_muscle(align):
        """Cette fonction renvoit un dict comprenant 3 choses :
        1. la séquence alignée
        2. la séquence sauvage
        3. la séquence avec les SNP
        """
        assert is_fasta(align), "Le fichier n'est pas un fichier fasta"

        def read_seq(sequence, index):
            """quick wrapper around AlignIO.read()"""
            return AlignIO.read(sequence, "fasta")[index]

        return {
            'exp': str(read_seq(align, 2).seq),  # la séquence expérimentale
            'snp': str(read_seq(align, 1).seq),  # le gene synthétique
            'wt': str(read_seq(align, 0).seq)  # la séquence de référence
        }

    def are_the_same(seq_list):
        """Détermine si les trois bases sont identiques. Renvoit `.` quand les trois bases sont
        identiques, `x` quand la base exp est la base `snp`, `X` quand la base exp est la base
        wt, et '-' quand la base exp ne correspond à aucune des possibilités

        """
        exp, snp, wt = seq_list

        if exp == '-':  # si la base est trimmée ou non alignée
            return '-'
        elif exp == snp and snp == wt:  # si les trois bases sont identiques
            return '.'
        elif exp == snp and snp != wt:  # si la base correspond au SNP
            return 'x'
        elif exp != snp and exp == wt:  # si la base correspond au WT
            return 'X'
        else:
            return 'N'  # si c'est encore autre chose ?

    def index_seq(sequence):
        """Détermine la position dans la séquence expérimentale"""

        position = 0
        position_list = []
        for index, base in enumerate(sequence):
            if base != '-': position += 1
            position_list.append(str(position))

        return position_list

    def polarity_snp(parsed_aln):
        """Give back the polarity of the SNP. Expect an output of parse_muscle."""

        consensus = []
        for i, base in enumerate(parsed_aln['exp']):
            exp = parsed_aln['exp'][i]
            snp = parsed_aln['snp'][i]
            wt = parsed_aln['wt'][i]

            consensus += are_the_same([exp, snp, wt])

        return consensus

    def parse_phd(align, rev=False):
        """Return a pandas DataFrame by parsing a phred output.

        :param align: a phd output.

        """

        def get_phd_seq_only(phd_file, reverse=False):
            """Return string from BEGIN_DNA to END_DNA. If reverse is True, give the
            sequence in inverse order."""
            ## inspiré de
            ## http://stackoverflow.com/questions/7559397/python-read-file-from-and-to-specific-lines-of-textq,
            ## réponse de EOL
            block = ""
            found = False
            with open(phd_file, 'r') as in_data:
                if not reverse:
                    for line in in_data:
                        if found:
                            if line.strip() == "END_DNA": break
                            else: block += line
                        else:
                            if line.strip() == "BEGIN_DNA":
                                found = True
                                block = ""
                else:
                    # if reverse flag is true, reverse sequence
                    for line in reversed(in_data.readlines()):
                        if found:
                            if line.strip() == "BEGIN_DNA": break
                            else: block += line
                        else:
                            if line.strip() == "END_DNA":
                                found = True
                                block = ""
            return block

        # pd.DataFrame.
        return pd.read_table(
            io.StringIO(u"%s" % get_phd_seq_only(align,
                                                 reverse=rev)),
            sep=" ",
            names=['base', 'qual', 'phase'])  # 'phase' = spectrogramme pos

    parsed_align = parse_muscle(input)

    clean_data = pd.DataFrame({
        'seqp': index_seq(parsed_align['exp']),
        'seqb': list(parsed_align['exp']),
        'refp': index_seq(parsed_align['wt']),
        'refb': list(parsed_align['wt']),
        'snpb': list(parsed_align['snp']),
        'cons': polarity_snp(parsed_align),
        'name': input
    })

    ## use pandas concat to concatenate two datasets along the axis 1, equivalent to cbind in R.
    ## ignore_index indicate that the two dataframes are not combined with their index.
    with_phd = pd.concat(
        [
            clean_data[clean_data.cons != '-'].reset_index(drop=True),
            # ne garde que les positions qui existent dans le phd file
            parse_phd(phd,
                      rev=reverse)
            # parse the phred output. reverse it if necessary to align to the reference.
        ],
        axis=1,
        ignore_index=True)
    # change les noms de colonnes.
    with_phd.columns = ['cons', 'name', 'refb', 'refp', 'seqb', 'seqp', 'snpb',
                        'base', 'qual', 'phase']

    if output != '-':
        with_phd.to_csv(output, index=False)
    else:
        click.echo(with_phd, file=output)

        ## tests

        # with open('pws1-1073bis.ab1.seq', 'rb') as pws:
        #     for record in SeqIO.parse(pws, 'fasta'):
        #         print record.seq.reverse_complement()
