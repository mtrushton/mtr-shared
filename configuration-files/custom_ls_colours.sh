# Customise the colours printed by the ls command. Below colours are assigned for each filetype. Note background colours can also be defined (this can be useful when dealing with dark background terminal schemes. 

# For directories: di=01;36;44 - This gives bold (01), cyan text (36) on a blue background (44)

#You can adjust the background color for directories by changing the last number:
#40 - Black background
#41 - Red background
#42 - Green background
#43 - Yellow background
#44 - Blue background
#45 - Magenta background
#46 - Cyan background
#47 - White background

#There are bright/light versions for all background colors, using the format 10X #where X is the normal color number:

#104 = Light blue (bright blue) background
#100 - Bright black (gray) background
#101 - Bright red background
#102 - Bright green background
#103 - Bright yellow background
#104 - Bright blue background
#105 - Bright magenta background
#106 - Bright cyan background
#107 - Bright white background

export LS_COLORS='di=01;30;104:fi=01;37:ex=01;91:ln=01;32:*.jpg=01;95:*.jpeg=01;95:*.png=01;95:*.gif=01;95:*.pdf=01;95:*.ps=01;95:*.txt=01;93:*.md=01;93:*.log=01;93:*.csv=01;93:*.json=01;93:*.xml=01;93:*.yml=01;93:*.yaml=01;93:*.html=01;93:*.css=01;93:*.js=01;93:*.zip=01;94:*.tar=01;94:*.gz=01;94:*.py=01;91:*.pro=01;91:*.f=01;91:*.sh=01;91:*.tex=01;93:*.fits=01;37'
