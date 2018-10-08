from pprint import pprint
from itertools import permutations
import sys

class Sequence:
    """
    ADT qui représente une séquence d'acide aminé
    """
    def __init__(self,sequence):
        self.seq= list(sequence)

    def __len__(self):
        return len(self.seq)

    def getItem(self,indice):
        return self.seq[indice]

    def setItem(self,indice,valeur):
        self.seq[indice] = valeur
        return self.seq


class Matrice:
    """
    ADT qui représente une matrice
    """

    def __init__(self,blosver):
        self.letters = {}
        self.mat = self.parse(blosver)

    def __getitem__(self,index):
        x,y = self.letters[index[0]],self.letters[index[1]]
        try:
            res = self.mat[x][y]
        except:
            res = self.mat[y][x]
        return res

    def parse(self,version):
        mat = []
        print(version)
        file = open(str(version), "r")
        if version[0] == "b":
            for i in range(6):
                line = file.readline()
        else:
            for i in range(9):
                line = file.readline()

        letters =file.readline().split()
        letters.pop()
        for i in range(len(letters)):
        	self.letters[letters[i]] = i

        while len(line.split()) > 0:
            line = file.readline()
            if line[0] == "*":
                break
            todel = list(line);del(todel[0]);del(todel[-1])
            line = "".join(todel)
            if len(line.split()) > 0:
                mat.append([int(i.strip()) for i in line.split()])
        return mat

class GlobalAligner:
    def __init__(self,blosPam,blos,V,W,seq1,seq2):
        self.S = blosPam
        self.V = V
        self.W = W
        self.blos = blos
        self.seq1 = seq1
        self.seq2 = seq2
        self.aligns = []

        self.align()

    def align(self,alignA=(),alignB=(),i=None,j=None):
        i = len(self.S)-1 if i == None else i
        j = len(self.S[0])-1 if j == None else j

        if (j > 0 or i > 0):
            if (i != 0 and self.S[i][j] == self.V[i][j]):
                self.align((self.seq1[i-1],)+alignA,("-",)+alignB,i-1,j)

            if (j != 0 and self.S[i][j] == self.W[i][j]):
                self.align(("-",)+alignA,(self.seq2[j-1],)+alignB,i,j-1)

            if (i > 0 and j > 0 and self.S[i][j] == self.S[i-1][j-1] + self.blos[self.seq1[i-1],self.seq2[j-1]]):
                self.align((self.seq1[i-1],)+alignA,(self.seq2[j-1],)+alignB,i-1,j-1)
        else:
            self.aligns.append([alignA,alignB])

def scoreMatrice(blosPam,seq1,seq2,e_gap,i_gap): #I = 1, E = 4

    V = Matrice(blosPam)
    V.matrix = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    W = Matrice(blosPam)
    W.matrix = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    S = Matrice(blosPam)
    S.matrix = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

    V.matrix[0][1] -= e_gap
    W.matrix[0][1] -= e_gap
    S.matrix[0][1] -= e_gap
    V.matrix[1][0] -= e_gap
    W.matrix[1][0] -= e_gap
    S.matrix[1][0] -= e_gap

    for i in range(2,len(V.matrix[0])):
    	V.matrix[0][i] = V.matrix[0][i-1] - i_gap
    	W.matrix[0][i] = W.matrix[0][i-1] - i_gap
    	S.matrix[0][i] = S.matrix[0][i-1] - i_gap

    for i in range(2,len(V.matrix)):
    	V.matrix[i][0] = V.matrix[i-1][0] - i_gap
    	W.matrix[i][0] = W.matrix[i-1][0] - i_gap
    	S.matrix[i][0] = S.matrix[i-1][0] - i_gap

    for i in range(1,len(V.matrix)):
        for j in range(1,len(V.matrix[0])):
            V.matrix[i][j] = max([S.matrix[i-1][j]-i_gap-e_gap, V.matrix[i-1][j]-e_gap])
            W.matrix[i][j] = max([S[i][j-1]-i_gap-e_gap, W.matrix[i][j-1]-e_gap])
            S.matrix[i][j] = max([S[i-1][j-1]+blos[seq1[i-1],seq2[j-1]], W.matrix[i][j],V.matrix[i][j]])
    return V,W,S

def findBestAligns(aligns,blosPam,e_gap,i_gap):
    """
        Fonction qui prend une liste d'alignements, la matrice BLOSUM et les gap initial and extended et va renvoyer
        les meilleurs alignements par rapport à ces paramètres
    """
    best_score = -1
    best_aligns = []

    for i in aligns:
        gap = False
        score = 0
        for j in range(len(i[0])):
            if not ("-" == i[0][j] or "-" == i[1][j]):
                score += blos[i[0][j],i[1][j]]
                gap = False
            else:
                if gap:
                    score -= e_gap
                else:
                    score -= i_gap
                    gap = True
        if score > best_score or best_score == -1:
            best_score = score
            best_aligns = [i]
        elif score == best_score:
            best_aligns += [i]
    return [best_score,best_aligns]


def seqParse(seqfile):
    """
    Renvoit l'ensemble de séquences d'un fichier
    """
    file = open(seqfile, "r")
    seqlist = []
    temp = ""
    for line in file:
        if line[0] != ">":

            temp += line.strip()
        else:
            if len(temp) > 0:
                seqlist.append(temp)
                temp = ""
    seqlist.append(Sequence(temp))
    return seqlist

def printAligns(score,aligns,scoreMat,length=70):
    print("Score:",score)
    res = [[] for i in aligns]

    for i in range(len(aligns)):
        alignment = str(i+1)
        alignment += "st" if i+1 == 1 else "nd" if i+1 == 2 else "rd" if i+1 == 3 else "th"
        alignment += " alignment"
        print(alignment+"-"*(length-len(alignment)))
        res[i] = [(aligns[i][0][k],aligns[i][1][k]) for k in range(len(aligns[i][0]))]

        for m in [res[i][k:k+length] if k+5 < len(res[i]) else i[k:len(res[i])-1] for k in range(0,len(res[i]),length)]:
            for j in [k[0] for k in m]:
                print(j,end="")
            print()
            for j in m:
                if (j[0] != "-" and j[1] != "-") and j[0] == j[1] :
                    print(":",end="")
                elif (j[0] != "-" and j[1] != "-") and scoreMat[j[0],j[1]] >= 0:
                    print(".",end="")
                else:
                    print(" ",end="")

            print()
            for j in [k[1] for k in m]:
                print(j,end="")
            print()
            print("-"*length)
        print()

if __name__ == '__main__':

    if len(sys.argv) == 3:
        filename = sys.argv[1]
        blosPam = sys.argv[2]
    else:
        sys.exit("Usage: " + sys.argv[0] + " [filename]")

    seqs=seqParse(filename)
    print("TEST")
    scoreMat = Matrice(blosPam)
    #scoreMat.initBlosPam(seqs)
    v,w,res = scoreMatrice(scoreMat,seqs[0],seqs[1],4,1)

    aligns = GlobalAligner(res,scoreMat,v,w,seqs[0],seqs[1])
    aligns.align()
    best_aligns = findBestAligns(aligns.aligns,scoreMat,4,1)
    printAligns(best_aligns[0],best_aligns[1],scoreMat)
