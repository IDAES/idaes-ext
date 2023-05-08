include ../Makefile.in

#ASL=/Users/john/src/idaes-ext-build/coinbrew/dist/include/coin-or/asl

ALL: functions.so

functions.o: functions.c
	$(CC) -c $(CFLAGS) functions.c -o functions.o -I$(ASL)

functions.so: functions.o
	$(LINK) $(LDFLAGS) functions.o -o functions.so

clean:
	rm -f *.o
	rm -f *.so
