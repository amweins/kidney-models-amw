	double precision abcis
	open (21,file='abcis.k')
	abcis=0.004
	do 10 i=1,26
        abcis=abcis + 0.001
  10    write(21,15) abcis
  15    format (e16.8)
	stop
	end

