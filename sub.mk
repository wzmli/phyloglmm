ms = makestuff
-include local.mk
-include $(ms)/os.mk

Sources += $(ms)

Makefile: $(ms)
$(ms):
	git submodule add https://github.com/dushoff/$@.git

$(ms)/%.mk: $(ms)
	git submodule init $(ms) 
	git submodule update $(ms) 
	touch $@
