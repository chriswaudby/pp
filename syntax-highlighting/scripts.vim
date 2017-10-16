if did_filetype()   " filetype already set..
  finish        " ..don't do these checks
endif
if search('#include <Avance\.incl>')
  setfiletype bruker
endif
