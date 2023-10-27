import pandas as pd

def GetOutput(genes, readcount):
    out = {}
    ENSG = list(genes.keys())
    out["ENSG"] = ENSG
    out["Gene"] = [genes[item].genename for item in ENSG]
    out["chrom"] = [genes[item].chrom for item in ENSG]
    out["start"] = [genes[item].start for item in ENSG]
    out["end"] = [genes[item].end for item in ENSG]
    out["strand"] = [genes[item].strand for item in ENSG]
    out["spliced"] = readcount['spliced'].squeeze().tolist()
    out["unspliced"] = readcount['unspliced'].squeeze().tolist()
    out["ambiguous"] = readcount['ambiguous'].squeeze().tolist()
    out["spanning"] = readcount['spanning'].squeeze().tolist()

    out = pd.DataFrame(out)

    return out