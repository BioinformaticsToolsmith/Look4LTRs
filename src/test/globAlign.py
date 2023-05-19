# A class that performs global alignment of two sequences
# It should have methods to return the length of the alignment, the actual alignments, the score, and the identity score

import numpy as np

# This is a class that performs global alignment of two sequences

class globAlign:
    def __init__(self, seq1, seq2, match=1, mismatch=-1, gap=-1):
        self.seq1 = seq1
        self.seq2 = seq2
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.score = 0
        self.identity = 0
        self.align1 = ""
        self.align2 = ""
        self.align()

    def align(self):
        # Initialize the matrices
        self.score = 0
        self.identity = 0
        self.align1 = ""
        self.align2 = ""
        self.scoreMatrix = np.zeros((len(self.seq1)+1, len(self.seq2)+1), dtype=int)
        self.directionMatrix = np.zeros((len(self.seq1)+1, len(self.seq2)+1), dtype=int)

        # Fill in the score matrix
        for i in range(1, len(self.seq1)+1):
            for j in range(1, len(self.seq2)+1):
                match = self.scoreMatrix[i-1, j-1] + self.match if self.seq1[i-1] == self.seq2[j-1] else self.scoreMatrix[i-1, j-1] + self.mismatch
                delete = self.scoreMatrix[i-1, j] + self.gap
                insert = self.scoreMatrix[i, j-1] + self.gap
                self.scoreMatrix[i, j] = max(0, match, delete, insert)
                if self.scoreMatrix[i, j] == 0:
                    self.directionMatrix[i, j] = 0
                elif self.scoreMatrix[i, j] == match:
                    self.directionMatrix[i, j] = 1
                elif self.scoreMatrix[i, j] == delete:
                    self.directionMatrix[i, j] = 2
                elif self.scoreMatrix[i, j] == insert:
                    self.directionMatrix[i, j] = 3

        # Traceback and compute the alignment
        i = len(self.seq1)
        j = len(self.seq2)
        score = 0
        identity = 0
        count = 0
        while i > 0 and j > 0:
            count += 1
            direction = self.directionMatrix[i, j]
            if direction == 0: # terminate
                break
            elif direction == 1: # match/mismatch
                self.align1 = self.seq1[i-1] + self.align1
                self.align2 = self.seq2[j-1] + self.align2
                if self.seq1[i-1] == self.seq2[j-1]:
                    identity += 1
                score += self.match if self.seq1[i-1] == self.seq2[j-1] else self.mismatch
                i -= 1
                j -= 1
            elif direction == 2: # delete
                self.align1 = self.seq1[i-1] + self.align1
                self.align2 = "-" + self.align2
                score += self.gap
                i -= 1
            elif direction == 3: # insert
                self.align1 = "-" + self.align1
                self.align2 = self.seq2[j-1] + self.align2
                score += self.gap
                j -= 1
        self.score = score
        self.identity = identity / count
    
    def getScore(self):
        return self.score

    def getIdentity(self):
        return self.identity

    def getAlign1(self):
        return self.align1

    def getAlign2(self):
        return self.align2

    def getLength(self):
        return len(self.align1)

    def printAlignment(self):
        print(self.align1)
        print(self.align2)

if __name__ == "__main__":
    seq1 = "AAGTAGGAAG"
    seq2 = "AAAAAAAAAA"
    align = globAlign(seq1, seq2)
    align.printAlignment()
    