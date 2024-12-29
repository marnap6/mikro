from Bio import SeqIO

prokka_file = "assembly.gbk"
bakta_file = "assembly.gbff"

def analyze_genes(file_path, file_format):
    gene_lengths = []
    polymerase_genes = []

    for record in SeqIO.parse(file_path, file_format):
        for feature in record.features:
            if feature.type == "CDS":
                if "translation" in feature.qualifiers:
                    gene_length = len(feature.qualifiers["translation"][0])
                    gene_lengths.append(gene_length)
                if "product" in feature.qualifiers:
                    product = feature.qualifiers["product"][0].lower()
                    gene_name = feature.qualifiers.get("gene", [None])[0]
                    if "rna polymerase" in product or "dna polymerase" in product:
                        if gene_name != None :
                            polymerase_genes.append(gene_name)
                        else:
                            polymerase_genes.append(product)
    gene_count = len(gene_lengths)
    avg_length = sum(gene_lengths) / gene_count if gene_count > 0 else 0
    return gene_count, avg_length, polymerase_genes

prokka_gene_count, prokka_avg_length, prokka_polymerases = analyze_genes(prokka_file, "genbank")
bakta_gene_count, bakta_avg_length, bakta_polymerases = analyze_genes(bakta_file, "genbank")

print("Porównanie wyników adnotacji:")
print(f"Prokka - Liczba genów: {prokka_gene_count}, Średnia długość genów: {prokka_avg_length:.2f}")
print(f"Bakta - Liczba genów: {bakta_gene_count}, Średnia długość genów: {bakta_avg_length:.2f}")

print("geny kodujace polimerazy RNA i DNA w Prokka:")
for gene in prokka_polymerases:
    print(f"- {gene}")
print(len(prokka_polymerases))

print("geny kodujace polimerazy RNA i DNA w Bakta:")
for gene in bakta_polymerases:
    print(f"- {gene}")
print(len(bakta_polymerases))
