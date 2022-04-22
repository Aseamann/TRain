# The following indices correspond to their respective CDR Loops, according to IMGT
c1 = [26,38]
c2 = [55,65]
c25 = [80,86]
c3 = [104,116]

class extract_CDRaa():
    # intialize extract_CDRaa
    def __init__(self, file):
        self.file = file

    # getSeqs returns entire AA sequences for each TCR gene family
    def getSeqs(self, lines, indices, a, flag):
        if flag == 0:
            start, end = indices[a], indices[a + 1]
            seq = lines[start:end]
        if flag == 1:
            start = indices[a]
            seq = lines[start:]
        i2 = [i for i, x in enumerate(seq[0]) if x == '|']
        gene = seq[0][i2[0] + 1:i2[1]]
        aa = seq[1:]
        aaSeq = '%s%s' % (aa[0][:-1], aa[1][:-1])
        return(gene, aaSeq)

    # getLoops pulls uses the CDR indices above to return a AA sequence
    def getLoops(self, seq_):
        cdr1 = seq_[c1[0]:c1[1]].replace('.', '')
        cdr2 = seq_[c2[0]:c2[1]].replace('.', '')
        cdr25 = seq_[c25[0]:c25[1]].replace('.', '')
        cdr3 = seq_[c3[0]:c3[1]].replace('.', '')
        whole = seq_.replace('.', '')
        return (cdr1, cdr2, cdr25, cdr3, whole)

    # storeSeqs produces a dictionary for TRAV and TRBV genes. Dictionary keys are
    # TR genes, the values are the AA sequences for each loop (c1, c2, c25, c3).
    def storeSeqs(self):
        imgtFile = open(self.file, 'r')
        cdrDict = {}
        lines = imgtFile.readlines()
        indices = [i for i, x in enumerate(lines) if x[0] == '>']
        for a in range(len(indices)):
            if a != len(indices) - 1:
                gene, aaSeq = self.getSeqs(lines, indices, a, 0)
                cd1, cd2, cd25, cd3, whole = self.getLoops(aaSeq)
                cdrDict[gene] = []
                cdrDict[gene].append(cd1)
                cdrDict[gene].append(cd2)
                cdrDict[gene].append(cd25)
                cdrDict[gene].append(cd3)
                cdrDict[gene].append(whole)
            else:
                gene, aaSeq = self.getSeqs(lines, indices, a, 1)
                cd1, cd2, cd25, cd3, whole = self.getLoops(aaSeq)
                cdrDict[gene] = []
                cdrDict[gene].append(cd1)
                cdrDict[gene].append(cd2)
                cdrDict[gene].append(cd25)
                cdrDict[gene].append(cd3)
                cdrDict[gene].append(whole)
        return(cdrDict)

# Call this function to retrieve AA sequence dictionaries
def cdr_loops():
    alp = extract_CDRaa('trav_alignment_imgt.txt')
    bet = extract_CDRaa('trbv_alignment_imgt.txt')
    alpha = alp.storeSeqs()
    beta = bet.storeSeqs()
    return(alpha, beta)

