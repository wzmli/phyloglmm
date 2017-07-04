ms = makestuff
-include $(ms)/os.mk

Sources += $(ms)

Makefile: $(ms)
$(ms):
	git submodule add https://github.com/dushoff/$@.git || mkdir $@

$(ms)/%.mk: 
	git submodule init $(ms) 
	git submodule update $(ms) 
	touch $@

