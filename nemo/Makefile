CFG=ORCA2_ICE_PISCES
ARCH=bisheng
MCFG=MY_ORCA2


.PHONY : all mkdir

all :
	./makenemo -m $(ARCH)  -j 32 -r $(CFG) -n '$(MCFG)' -d 'OCE' del_key 'key_si3 key_top key_xios key_qco key_linssh' add_key 'MTHREAD COUPLE'



clean :
	./makenemo -m $(ARCH) -r $(CFG) -n '$(MCFG)' -j 4 clean

