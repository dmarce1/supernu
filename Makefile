all:	
	cd src && make && cp supernu ..
clean: 
	cd src && make clean && rm ../supernu

