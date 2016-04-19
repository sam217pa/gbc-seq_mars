  #!/usr/bin/env python

  from Bio import SeqIO
  from glob import glob
  from os.path import basename

  def sw_or_ws(mutant_name):
      """
      Determine si le mutant est SW ou WS
      """
      if 'sw' in mutant_name:
          return 'sw'
      else:
          return 'ws'

  # en tete de colonne
  print "id name mutant"

  def print_seq_id(dna_seq):
      with open(dna_seq, "rb") as spectro:
          sequence = list(SeqIO.parse(spectro, "abi"))
          print sequence[0].id + " " + sequence[0].name + " " + sw_or_ws(sequence[0].name)

  for file in glob("data/sw/spectro/*.ab1"):
      if basename(file) != "psw76-1073bis.ab1": # exclut la sequence qui pose probleme
          print_seq_id(file)

  for file in glob("data/ws/spectro/*.ab1"):
      print_seq_id(file)
