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

def readFromFastaOrdered(filePath, removeDashes = False):
    sequences = []
    currentSequence = None

    with open(filePath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):                    
                tag = line[1:]
                currentSequence = Sequence(tag, "")
                sequences.append(currentSequence)
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

#reads match columns only
def readFromStockholm(filePath, includeInsertions = False):
    sequences = {}
    
    with open(filePath, 'r') as stockFile:
        for line in stockFile:
            line = line.strip()
            if line == "//":
                break
            elif line == "" or line[0] == "#":
                pass
            else:  
                key, seq = line.split()
                if key not in sequences:
                    sequences[key] = Sequence(key, "")
                    
                for c in seq:
                    #if includeInsertions or not (c == '.' or c in string.ascii_lowercase):
                    if includeInsertions or (c == c.upper() and c != '.'):
                        sequences[key].seq = sequences[key].seq + c    
    return sequences

def writeFasta(alignment, filePath, taxa = None, append = False):
        with open(filePath, 'a' if append else 'w') as textFile:
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

def cleanGapColumns(filePath, cleanFile = None):
    align = readFromFasta(filePath, False)
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
        cleanFile = filePath
        
    writeFasta(align, cleanFile)
    
def convertRnaToDna(filePath, destFile = None):
    align = readFromFasta(filePath, False)
    for taxon in align:
        align[taxon].seq = align[taxon].seq.replace('U', 'T')
    if destFile is None:
        destFile = filePath
    writeFasta(align, destFile)

def inferDataType(filePath):
    sequences = readFromFasta(filePath, removeDashes=True)
    acg, t, u, total = 0, 0, 0, 0
    for taxon in sequences:
        letters = sequences[taxon].seq.upper()
        for letter in letters:
            total = total + 1
            
            if letter in ('A', 'C', 'G', 'N'):
                acg = acg + 1
            elif letter == 'T':
                t = t + 1
            elif letter == 'U':
                u = u + 1
    
    if u == 0 and (acg + t)/total > 0.9:
        print("Found {}% ACGT-N, assuming DNA..".format(int(100*(acg + t)/total)))
        dataType = "dna"
    elif t == 0 and (acg + u)/total > 0.9:
        print("Found {}% ACGU-N, assuming RNA..".format(int(100*(acg + u)/total)))
        dataType = "rna"
    else:
        print("Assuming protein..")
        dataType = "protein"
          
    return dataType

def readSequenceLengthFromFasta(filePath):
    with open(filePath) as f:
        length = 0
        readSequence = False
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if readSequence:
                    return length
                readSequence = True
            else:
                length = length + len(line)
    if readSequence:
        return length

def countGaps(alignFile):
    counts = []
    currentSequence = ""

    with open(alignFile) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'): 
                
                if currentSequence is not None:
                    if len(counts) == 0:
                        counts = [0] * len(currentSequence)
                    for i in range(len(counts)):
                        if currentSequence[i] == '-':
                            counts[i] = counts[i] + 1
                                             
                currentSequence = ""
            else:
                currentSequence = currentSequence + line
        if currentSequence is not None:
            if len(counts) == 0:
                counts = [0] * len(currentSequence)
            for i in range(len(counts)):
                if currentSequence[i] == '-':
                    counts[i] = counts[i] + 1
    
    return counts
