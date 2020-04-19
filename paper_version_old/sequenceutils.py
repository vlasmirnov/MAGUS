'''
Created on Sep 22, 2018

@author: Vlad
'''


class Sequence:
    def __init__(self, tag, seq):
        self.tag = tag
        self.seq = seq
    
def readFromFasta(filePath, removeDashes = False):
    sequences = {}
    currentSequence = None

    with open(filePath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):                    
                tag = line[1:]
                currentSequence = Sequence(tag, "")
                sequences[tag] = currentSequence
            else :
                if(removeDashes):
                    line = line.replace("-", "")
                currentSequence.seq = currentSequence.seq + line

    print("Read " + str(len(sequences)) + " sequences from " + filePath + " ..")
                                
    return sequences


def readFromPhylip(filePath, removeDashes = False):
    sequences = {}    

    with open(filePath) as f:
        firstLine = f.readline().strip()
        
        for line in f:
            #print(line)
            tokens = line.split()
            if len(tokens) == 2:
                tag = tokens[0]
                seq = tokens[1]
                
                if(removeDashes):
                    seq = seq.replace("-", "")
                   
                if tag in sequences:
                    sequences[tag].seq = sequences[tag].seq + seq
                else:
                    sequences[tag] = Sequence(tag, seq)
                
    
    print("Read " + str(len(sequences)) + " sequences from " + filePath + " ..")
                                
    return sequences


def writeFasta(alignment, filePath, taxa = None):
        with open(filePath, 'w') as textFile:
            if taxa is not None:
                for tag in taxa:
                    if tag in alignment:
                        textFile.write('>' + tag + '\n' + alignment[tag].seq + '\n')
            else:
                for tag in alignment:
                    textFile.write('>' + tag + '\n' + alignment[tag].seq + '\n')
      
                    
def writePhylip(alignment, filePath, taxa = None):
    maxChars = 0
    lines = []
    for tag in alignment:
        if taxa is None or tag in taxa:
            lines.append("{} {}\n".format(tag, alignment[tag].seq))
            maxChars = max(maxChars, len(alignment[tag].seq))
    
    with open(filePath, 'w') as textFile:
        textFile.write("{} {}\n".format(len(lines), maxChars))
        for line in lines:
            textFile.write(line)

def cleanGapColumns(alignFile, cleanFile = None):
    align = readFromFasta(alignFile, False)
    values = list(align.values())
    keepCols = []
    for i in range(len(values[0].seq)):
        for j in range(len(values)):
            if values[j].seq[i] != '-':
                keepCols.append(i)
                break
            
    print("Removing gap columns.. Kept {} out of {}..".format(len(keepCols), len(values[0].seq)))
    for s in values:
        s.seq = ''.join(s.seq[idx] for idx in keepCols)
    
    if cleanFile is None:
        cleanFile = alignFile
        
    writeFasta(align, cleanFile)
    
def convertDnaToRna(srcFile, destFile = None):
    align = readFromFasta(srcFile, False)
    for taxon in align:
        align[taxon].seq = align[taxon].seq.replace('U', 'T')
    if destFile is None:
        destFile = srcFile
    writeFasta(align, destFile)
            
