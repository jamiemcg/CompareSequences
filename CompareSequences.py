#!/usr/bin/env python3

# CompareSequences.py
# Jamie McGowan, 2019 <jamie.mcgowan@mu.ie>
# Usage: python CompareSequences.py one.fasta two.fasta


import sys
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sequences_a = []
sequences_b = []

nonpolar = ["G", "A", "V", "L", "I", "M", "F", "W", "P"]
polar = ["S", "T", "C", "Y", "N", "Q"]
negative = ["D", "E"]
positive = ["K", "R", "H"]

amino_acids = nonpolar + polar + negative + positive

for record in SeqIO.parse(sys.argv[1], "fasta"):
	sequences_a.append(record)

for record in SeqIO.parse(sys.argv[2], "fasta"):
	sequences_b.append(record)

# Bar chart showing number of sequences in each file
number_of_sequences = pd.DataFrame([[sys.argv[1], len(sequences_a)], [sys.argv[2], len(sequences_b)]], columns = ["Filename", "Number of Sequences"])

ax = sns.barplot(x = "Filename", y = "Number of Sequences", data = number_of_sequences)
plt.savefig("plotNumberOfSequences.pdf")
plt.clf()
# plt.show()

# Box plot showing the number of each amino acid per protein
for aa in amino_acids:
	aa_three_letter = SeqUtils.seq3(aa)
	
	amino_acids_count = []

	for protein in sequences_a:
		amino_acids_count.append([sys.argv[1], str(protein.seq).count(aa)])

	for protein in sequences_b:
		amino_acids_count.append([sys.argv[2], str(protein.seq).count(aa)])

	amino_acids_count = pd.DataFrame(amino_acids_count, columns = ["Filename", "Number of " + aa_three_letter + " residues"])

	ax = sns.boxplot(x = "Filename", y = "Number of " + aa_three_letter + " residues", data = amino_acids_count)
	plt.savefig("plotNumberOf" + aa_three_letter + "Residues.pdf")
	plt.clf()
	# plt.show()

# Box plot showing the percentage of each amino acid per protein
for aa in amino_acids:
	aa_three_letter = SeqUtils.seq3(aa)
	
	amino_acids_percentage = []

	for protein in sequences_a:
		amino_acids_percentage.append([sys.argv[1], (100 * (1.0 * str(protein.seq).count(aa)) / len(str(protein.seq)))])

	for protein in sequences_b:
		amino_acids_percentage.append([sys.argv[2], (100 * (1.0 * str(protein.seq).count(aa)) / len(str(protein.seq)))])

	amino_acids_percentage = pd.DataFrame(amino_acids_percentage, columns = ["Filename", "Percentage of " + aa_three_letter + " residues"])

	ax = sns.boxplot(x = "Filename", y = "Percentage of " + aa_three_letter + " residues", data = amino_acids_percentage)
	plt.savefig("plotPercentageOf" + aa_three_letter + "Residues.pdf")
	plt.clf()
	# plt.show()

# Box plot showing distribution of protein lengths
protein_lengths = []

for protein in sequences_a:
	protein_lengths.append([sys.argv[1], len(str(protein.seq))])

for protein in sequences_b:
	protein_lengths.append([sys.argv[2], len(str(protein.seq))])

protein_lengths = pd.DataFrame(protein_lengths, columns = ["Filename", "Protein lengths"])

ax = sns.boxplot(x = "Filename", y = "Protein lengths", data = protein_lengths)
plt.savefig("plotProteinLengths.pdf")
plt.clf()
# plt.show()

# Box plot showing distribution of protein molecular weights
molecular_weight = []

for protein in sequences_a:
	molecular_weight.append([sys.argv[1], SeqUtils.molecular_weight(str(protein.seq).replace("X", ""), "protein")])

for protein in sequences_b:
	molecular_weight.append([sys.argv[2], SeqUtils.molecular_weight(str(protein.seq).replace("X", ""), "protein")])

molecular_weight = pd.DataFrame(molecular_weight, columns = ["Filename", "Molecular Weight"])

ax = sns.boxplot(x = "Filename", y = "Molecular Weight", data = molecular_weight)
plt.savefig("plotMolecularWeights.pdf")
plt.clf()
# plt.show()

gravy_index = []
aromaticity = []
instability_index = []
# flexibility = []
isoelectric_point = []
secondary_structure_fraction = []

for protein in sequences_a:
	analysed_seq = ProteinAnalysis(str(protein.seq).replace("X", ""))

	gravy_index.append([sys.argv[1], analysed_seq.gravy()])
	aromaticity.append([sys.argv[1], analysed_seq.aromaticity()])
	instability_index.append([sys.argv[1], analysed_seq.instability_index()])
	# flexibility.append([sys.argv[1], analysed_seq.flexibility()])
	isoelectric_point.append([sys.argv[1], analysed_seq.isoelectric_point()])
	secondary_structure_fraction.append([sys.argv[1], analysed_seq.secondary_structure_fraction()])

for protein in sequences_b:
	analysed_seq = ProteinAnalysis(str(protein.seq).replace("X", ""))

	gravy_index.append([sys.argv[2], analysed_seq.gravy()])
	aromaticity.append([sys.argv[2], analysed_seq.aromaticity()])
	instability_index.append([sys.argv[2], analysed_seq.instability_index()])
	# flexibility.append([sys.argv[2], analysed_seq.flexibility()])
	isoelectric_point.append([sys.argv[2], analysed_seq.isoelectric_point()])
	secondary_structure_fraction.append([sys.argv[2], analysed_seq.secondary_structure_fraction()])

# Box plot showing gravy indexes
gravy_index = pd.DataFrame(gravy_index, columns = ["Filename", "Gravy Index"])

ax = sns.boxplot(x = "Filename", y = "Gravy Index", data = gravy_index)
plt.savefig("plotGravyIndex.pdf")
# plt.show()
plt.clf()

aromaticity = pd.DataFrame(aromaticity, columns = ["Filename", "Aromaticity"])

ax = sns.boxplot(x = "Filename", y = "Aromaticity", data = aromaticity)
plt.savefig("plotAromaticity.pdf")
# plt.show()
plt.clf()

instability_index = pd.DataFrame(instability_index, columns = ["Filename", "Instability Index"])

ax = sns.boxplot(x = "Filename", y = "Instability Index", data = instability_index)
plt.savefig("plotInstabilityIndex.pdf")
# plt.show()
plt.clf()

isoelectric_point = pd.DataFrame(isoelectric_point, columns = ["Filename", "Isoelectric Point"])

ax = sns.boxplot(x = "Filename", y = "Isoelectric Point", data = isoelectric_point)
plt.savefig("plotIsoelectricPoint.pdf")
# plt.show()
plt.clf()

test_secondary_structure_fraction = []

# Reformat data for secondary structure fractions to generate a data frame
for i in range(0, len(secondary_structure_fraction)):
	filename = secondary_structure_fraction[i][0]
	helix = secondary_structure_fraction[i][1][0]
	turn = secondary_structure_fraction[i][1][1]
	sheet = secondary_structure_fraction[i][1][2]

	# secondary_structure_fraction[i] = [filename, helix, turn, sheet]
	
	test_secondary_structure_fraction.append([filename, "Helix", helix])
	test_secondary_structure_fraction.append([filename, "Turn", turn])
	test_secondary_structure_fraction.append([filename, "Sheet", sheet])

# secondary_structure_fraction = pd.DataFrame(secondary_structure_fraction, columns = ["Filename", "Helix", "Turn", "Sheet"])

test_secondary_structure_fraction = pd.DataFrame(test_secondary_structure_fraction, columns = ["Filename", "Structure", "Fraction"])

ax = sns.violinplot(x = "Filename", y = "Fraction", hue = "Structure", data = test_secondary_structure_fraction)
plt.savefig("plotSecondaryStructureFraction.pdf")
# plt.show()
plt.clf()

ax = sns.violinplot(x = "Structure", y = "Fraction", hue = "Filename", data = test_secondary_structure_fraction)
plt.savefig("plotSecondaryStructureFractionAlt.pdf")
# plt.show()
plt.clf()
