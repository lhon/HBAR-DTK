unitigger = bogart
utgErrorRate = 0.015 
utgErrorLimit = 4.5

cnsErrorRate = 0.25
cgwErrorRate = 0.25
ovlErrorRate = 0.015

frgMinLen = 1000
ovlMinLen = 40
  
merSize=14

merylMemory = 1311
merylThreads = 8

ovlStoreMemory = 1311

# grid info
useGrid = 0
scriptOnGrid = 0
frgCorrOnGrid = 0
ovlCorrOnGrid = 0

sge = -S /bin/bash -sync y -V -q smrtpipe-fast
sgeScript = -pe smp 6
sgeConsensus = -pe smp 1
sgeOverlap = -pe smp 1
sgeFragmentCorrection = -pe smp 2
sgeOverlapCorrection = -pe smp 1

#ovlHashBits = 22
#ovlHashBlockLength = 46871347
#ovlRefBlockSize =  537

ovlHashBits = 24
ovlThreads = 3
ovlHashBlockLength = 20000000
ovlRefBlockSize =  50000000



#ovlThreads = 1
ovlConcurrency = 6
frgCorrThreads = 2
frgCorrBatchSize = 100000
ovlCorrBatchSize = 100000

sgeName = 056518
