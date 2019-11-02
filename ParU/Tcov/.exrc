if &cp | set nocp | endif
nnoremap  
let s:cpo_save=&cpo
set cpo&vim
nnoremap <NL> <NL>
nnoremap  
nnoremap  
nnoremap <silent>  :CtrlP
nnoremap <silent> w :CCTreeWindowToggle
nnoremap <silent> y :CCTreeWindowSaveCopy
nmap vd :vert scs find d =expand("<cword>")
nmap vi :vert scs find i ^=expand("<cfile>")$	
nmap vf :vert scs find f =expand("<cfile>")	
nmap ve :vert scs find e =expand("<cword>")
nmap vt :vert scs find t =expand("<cword>")
nmap vc :vert scs find c =expand("<cword>")
nmap vg :vert scs find g =expand("<cword>")
nmap vs :vert scs find s =expand("<cword>")
nmap hd :scs find d =expand("<cword>")	
nmap hi :scs find i ^=expand("<cfile>")$	
nmap hf :scs find f =expand("<cfile>")	
nmap he :scs find e =expand("<cword>")	
nmap ht :scs find t =expand("<cword>")	
nmap hc :scs find c =expand("<cword>")	
nmap hg :scs find g =expand("<cword>")	
nmap hs :scs find s =expand("<cword>")	
nmap d :cs find d =expand("<cword>")	
nmap i :cs find i ^=expand("<cfile>")$
nmap f :cs find f =expand("<cfile>")	
nmap e :cs find e =expand("<cword>")	
nmap t :cs find t =expand("<cword>")	
nmap c :cs find c =expand("<cword>")	
nmap g :cs find g =expand("<cword>")	
nmap s :cs find s =expand("<cword>")	
nnoremap \  :nohlsearch
nmap gx <Plug>NetrwBrowseX
nnoremap j gj
nnoremap k gk
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#NetrwBrowseX(expand("<cWORD>"),0)
nnoremap <SNR>62_: :=v:count ? v:count : ''
xnoremap <silent> <Plug>(Limelight) :Limelight
nnoremap <silent> <Plug>(Limelight) :set opfunc=limelight#operatorg@
nnoremap <silent> <Plug>(vimshell_create) :VimShellCreate
nnoremap <silent> <Plug>(vimshell_switch) :VimShell
nmap <Nul><Nul>d :vert scs find d =expand("<cword>")
nmap <Nul><Nul>i :vert scs find i ^=expand("<cfile>")$	
nmap <Nul><Nul>f :vert scs find f =expand("<cfile>")	
nmap <Nul><Nul>e :vert scs find e =expand("<cword>")
nmap <Nul><Nul>t :vert scs find t =expand("<cword>")
nmap <Nul><Nul>c :vert scs find c =expand("<cword>")
nmap <Nul><Nul>g :vert scs find g =expand("<cword>")
nmap <Nul><Nul>s :vert scs find s =expand("<cword>")
nmap <Nul>d :scs find d =expand("<cword>")	
nmap <Nul>i :scs find i ^=expand("<cfile>")$	
nmap <Nul>f :scs find f =expand("<cfile>")	
nmap <Nul>e :scs find e =expand("<cword>")	
nmap <Nul>t :scs find t =expand("<cword>")	
nmap <Nul>c :scs find c =expand("<cword>")	
nmap <Nul>g :scs find g =expand("<cword>")	
nmap <Nul>s :scs find s =expand("<cword>")	
map <F2> :NERDTreeToggle 
let &cpo=s:cpo_save
unlet s:cpo_save
set autoindent
set backspace=indent,eol,start
set backup
set backupdir=~/.vim/.backup//
set cscopeprg=/usr/bin/cscope
set cscopetag
set cscopeverbose
set directory=~/.vim/.swp//
set errorbells
set expandtab
set fileencodings=ucs-bom,utf-8,latin1
set guicursor=n-v-c:block,o:hor50,i-ci:hor15,r-cr:hor30,sm:block,a:blinkon0
set helplang=en
set hlsearch
set incsearch
set lazyredraw
set mouse=a
set ruler
set runtimepath=~/.vim,~/.vim/bundle/vimshell.vim,~/.vim/bundle/vimproc.vim,~/.vim/bundle/cscopemaps.vim,~/.vim/bundle/matlab,~/.vim/bundle/goyo.vim,~/.vim/bundle/limelight.vim,~/.vim/bundle/Vundle.vim,~/.vim/bundle/nerdtree,~/.vim/bundle/vim-airline,~/.vim/bundle/vim-fugitive,~/.vim/bundle/vimux,~/.vim/bundle/vim-mac-classic-theme,~/.vim/bundle/vim-colors-solarized,~/.vim/bundle/ctrlp.vim,/usr/share/vim/vimfiles,/usr/share/vim/vim74,/usr/share/vim/vimfiles/after,~/.vim/after,~/.vim/bundle/Vundle.vim,~/.vim/bundle/vimshell.vim/after,~/.vim/bundle/vimproc.vim/after,~/.vim/bundle/cscopemaps.vim/after,~/.vim/bundle/matlab/after,~/.vim/bundle/goyo.vim/after,~/.vim/bundle/limelight.vim/after,~/.vim/bundle/Vundle.vim/after,~/.vim/bundle/nerdtree/after,~/.vim/bundle/vim-airline/after,~/.vim/bundle/vim-fugitive/after,~/.vim/bundle/vimux/after,~/.vim/bundle/vim-mac-classic-theme/after,~/.vim/bundle/vim-colors-solarized/after,~/.vim/bundle/ctrlp.vim/after,~/.vim/eclim,~/.vim/eclim/after
set shell=/bin/bash\ -l
set shiftwidth=4
set showcmd
set showmatch
set smartindent
set softtabstop=4
set tabstop=4
set tags=./tags;/users/aznaveh
set visualbell
set wildmenu
" vim: set ft=vim :
