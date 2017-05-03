" Vim syntax file
" Language: bruker pulse programs
" Maintainer: Chris Waudby
" Latest Revision: 3 March 2017

if exists("b:current_syntax")
    finish
endif

syn keyword basicLanguageKeywords center larger lalign ralign go BLKGRAD UNBLKGRAD zd ze mc prosol fq
syn keyword brukerCommand define if else goto print define  
syn match brukerPower "pl[0-9]+"
syn match brukerChannel /:f[1-4]/ms=s+1


syn region brukerExpr start='"' end='"'
syn region brukerExpr start='<' end='>'

syn match brukerComment ";.*$"
syn region brukerComment start='/\*' end='\*/'

syn match brukerPrecomp "#.*$"

let b:current_syntax = "bruker"

hi def link brukerComment       Comment
hi def link brukerPrecomp       PreProc
hi def link brukerExpr          Constant
hi def link brukerCommand       Statement
hi def link brukerChannel       Type
hi def link brukerPower         PreProc
hi def link basicLanguageKeywords   Identifier
