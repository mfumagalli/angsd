B1;2cCC=gcc
CXX=g++

FLAGS=-O3

all: angsd misc

.PHONY: misc


misc:
	make -C misc/



bfgs.o: bfgs.cpp bfgs.h
	$(CXX) $(FLAGS) -c bfgs.cpp -Wno-write-strings

analysisFunction.o: analysisFunction.cpp analysisFunction.h shared.o
	$(CXX) $(FLAGS) -c analysisFunction.cpp

analysisAsso.o: analysisFunction.o analysisAsso.cpp analysisAsso.h shared.o
	$(CXX) $(FLAGS) -c analysisAsso.cpp

analysisMajorMinor.o: analysisFunction.o analysisMajorMinor.cpp analysisMajorMinor.h shared.o
	$(CXX) $(FLAGS) -c analysisMajorMinor.cpp


analysisFasta.o: analysisFunction.o analysisFasta.cpp analysisFasta.h shared.o
	$(CXX) $(FLAGS) -c analysisFasta.cpp

analysisAbbababa.o: analysisFunction.o analysisAbbababa.cpp analysisAbbababa.h shared.o
	$(CXX) $(FLAGS) -c analysisAbbababa.cpp


analysisEstLikes.o: analysisFunction.o analysisEstLikes.cpp analysisEstLikes.h shared.o
	$(CXX) $(FLAGS) -c analysisEstLikes.cpp


analysisKeepList.o: analysisKeepList.h analysisKeepList.cpp analysisFunction.o shared.o keepToBin.cpp bams.h
	$(CXX) $(FLAGS) -c analysisKeepList.cpp


analysisMaf.inbreed.o: analysisFunction.o analysisMaf.inbreed.cpp analysisMaf.h shared.o analysisFunction.o bfgs.o
	$(CXX) $(FLAGS) -c analysisMaf.inbreed.cpp -c

analysisMaf.o: analysisFunction.o analysisMaf.cpp analysisMaf.h shared.o analysisFunction.o bfgs.o
	$(CXX) $(FLAGS) -c analysisMaf.cpp

analysisCount.o: analysisCount.cpp analysisCount.h analysisFunction.o general.h kstring.o
	$(CXX) $(FLAGS) -c analysisCount.cpp

analysisAncError.o: analysisAncError.cpp analysisAncError.h analysisFunction.o
	$(CXX) $(FLAGS) -c analysisAncError.cpp


kprobaln.o: kprobaln.c  kprobaln.h
	$(CC) $(FLAGS) -c kprobaln.c

baq_adjustMapQ.o: baq_adjustMapQ.cpp baq_adjustMapQ.h kprobaln.o
	$(CC) $(FLAGS) -c baq_adjustMapQ.cpp




general.o: general.cpp general.h shared.o analysisCovar.cpp analysisEstError.o analysisCallGenotypes.o  angsd_realSFS.cpp analysisEstLikes.o analysisMajorMinor.o  analysisCount.o analysisAsso.o analysisFunction.o analysisKeepList.cpp  analysisKeepList.o thorfinn.o snpStat.o analysisHWE.o analysisAncError.o getFasta.o snptools.cpp hetplas.cpp writePlink.cpp analysisAbbababa.o analysisFasta.o
	$(CXX) $(FLAGS) -c general.cpp


reader_sim.o: shared.o reader_sim.cpp reader_sim.h
	$(CXX) $(FLAGS) -c reader_sim.cpp
 
pileup.o: pileup.cpp pileup.h 
	$(CXX) $(FLAGS) -c pileup.cpp



analysisEstError.o: analysisEstError.cpp analysisEstError.h analysisFunction.h analysisMaf.h shared.o
	$(CXX) $(FLAGS) -c analysisEstError.cpp


shared.o: shared.cpp shared.h 
	$(CXX) $(FLAGS) -c shared.cpp 
printRes.o: shared.o printRes.cpp printRes.h
	$(CXX) $(FLAGS) -c printRes.cpp

soapReader.o: soapReader.cpp soapReader.h
	$(CXX) $(FLAGS) -c soapReader.cpp

getFasta.o: getFasta.cpp getFasta.h
	$(CXX) $(FLAGS) -c getFasta.cpp

beagleReader.o: beagleReader.cpp beagleReader.h
	$(CXX) $(FLAGS) -c beagleReader.cpp


soapMaster.o: soapMaster.cpp soapMaster.h shared.h soapReader.o analysisFunction.o
	$(CXX) $(FLAGS) -c soapMaster.cpp

tglfs.o: tglfs.cpp tglfs.h 
	$(CXX) $(FLAGS) -c tglfs.cpp 

multiReader.o: multiReader.cpp multiReader.h soapMaster.o soapReader.o pileups.o pileup.o beagleReader.o tglfs.o reader_sim.o
	$(CXX) $(FLAGS) -c multiReader.cpp


pileups.o: pileups.h pileups.cpp pileup.o
	$(CXX) $(FLAGS) -c pileups.cpp

razf.o: razf.c razf.h
	$(CC)  $(FLAGS) razf.c -c -D_USE_KNETFILE

knetfile.o: knetfile.c knetfile.h
	$(CC)  $(FLAGS) knetfile.c -c


faidx.o: faidx.c faidx.h razf.o
	$(CC)  $(FLAGS) faidx.c -c -D_USE_KNETFILE

bgzf.o:	bgzf.c bgzf.h knetfile.o
	$(CC) $(FLAGS) bgzf.c -c -D_USE_KNETFILE

kstring.o: kstring.h kstring.c
	$(CC) $(FLAGS) kstring.c -c 

bams.o : bams.cpp bams.h baq_adjustMapQ.o
	$(CXX) bams.cpp $(FLAGS) -c -Wunused

indexer.o : bams.h indexer.h indexer.cpp
	$(CXX) indexer.cpp $(FLAGS) -c

mUpPile.o: mUpPile.cpp mUpPile.h bams.o threadpool2.cpp faidx.o bambi_interface.h
	$(CXX) mUpPile.cpp $(FLAGS) -c

parseArgs_bambi.o: parseArgs_bambi.cpp parseArgs_bambi.h
	$(CXX) parseArgs_bambi.cpp $(FLAGS) -c

bam_likes.o: bam_likes.cpp bam_likes.h
	$(CXX) bam_likes.cpp $(FLAGS) -c

soap_likes.o: soap_likes.cpp soap_likes.h
	$(CXX) soap_likes.cpp $(FLAGS) -c # -fomit-frame-pointer  -ffast-math -funroll-loops -mmmx -msse -msse2 -msse3 -fmessage-length=0 -maccumulate-outgoing-args -MMD -MP 

thorfinn.o: general.h thorfinn.cpp thorfinn.h
	$(CXX) $(FLAGS) -c thorfinn.cpp  

gatk_likes.o: gatk_likes.cpp gatk_likes.h
	$(CXX) gatk_likes.cpp $(FLAGS) -c

bammer_main.o : parseArgs_bambi.o bammer_main.cpp mUpPile.o bams.o kstring.o indexer.o threadpool2.cpp faidx.o
	$(CXX) bammer_main.cpp $(FLAGS) -c


analysisHWE.o: analysisFunction.h shared.h analysisHWE.cpp analysisHWE.h
	$(CXX) analysisHWE.cpp $(FLAGS) -c

analysisCallGenotypes.o: analysisFunction.h shared.h analysisCallGenotypes.cpp analysisCallGenotypes.h
	$(CXX) analysisCallGenotypes.cpp $(FLAGS) -c


snpStat.o: analysisFunction.h shared.h snpStat.cpp snpStat.h analysisHWE.h
	$(CXX) snpStat.cpp $(FLAGS) -c


angsd: angsd.cpp shared.o bfgs.o multiReader.o printRes.o general.o soap_likes.o gatk_likes.o bam_likes.o bammer_main.o analysisFunction.o bgzf.o analysisMaf.o
	$(CXX) $(FLAGS) angsd.cpp *.o  -lpthread -lz -o angsd

inbreed: angsd.cpp shared.o bfgs.o multiReader.o printRes.o general.o soap_likes.o gatk_likes.o bam_likes.o bammer_main.o analysisFunction.o bgzf.o analysisMaf.inbreed.o
	$(CXX) $(FLAGS) angsd.cpp *.o  -lpthread -lz -o angsd_inbreed

test:
	@echo "angsd: test scripts not implemented yet."

clean:
	@rm  -f *.o angsd *~ angsd_inbreed
