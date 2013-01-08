# Test with an old compiler (Apple gcc 4.2)
CC=/usr/bin/g++ -Wall
# Test with a mondern compiler (Apple clang 4.1)
#CC=clang++ -Wall

all: Rolke FeldmanCousins test_TMath

SpecFuncCephes.o: SpecFuncCephes.h SpecFuncCephes.cxx
	$(CC) -c SpecFuncCephes.cxx

TMath.o: TMath.h TMath.cxx SpecFuncCephes.o
	$(CC) -c TMath.cxx

TRolke.o: TRolke.h TRolke.cxx TMath.o
	$(CC) -c TRolke.cxx

Rolke.o: Rolke.C TRolke.o
	$(CC) -c Rolke.C

Rolke: Rolke.o TRolke.o
	$(CC) Rolke.o TRolke.o TMath.o SpecFuncCephes.o -lm -o Rolke

TFeldmanCousins.o: TFeldmanCousins.h TFeldmanCousins.cxx TMath.o
	$(CC) -c TFeldmanCousins.cxx

FeldmanCousins.o: FeldmanCousins.C TFeldmanCousins.o
	$(CC) -c FeldmanCousins.C

FeldmanCousins: FeldmanCousins.o TFeldmanCousins.o
	$(CC) FeldmanCousins.o TFeldmanCousins.o TMath.o SpecFuncCephes.o -lm -o FeldmanCousins

test_TMath.o: test_TMath.C TMath.o
	$(CC) -c test_TMath.C

test_TMath: test_TMath.o SpecFuncCephes.o
	$(CC) test_TMath.o TMath.o SpecFuncCephes.o -lm -o test_TMath

clean:
	- rm *.o Rolke FeldmanCousins test_TMath
