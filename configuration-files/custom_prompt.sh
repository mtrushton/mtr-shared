# Custom prompts in Bash. Update the .bashrc file with one of the following or include source custom_prompt.sh

# Kali linux prompt running in Bash
#PS1=$'\\[\e[1;31m\\]\u250c\u2500\u2500(\\[\e[1;33m\\]\\u@\\h\\[\e[1;31m\\])\u2500[\\[\e[0;37m\\]\\w\\[\e[1;31m\\]]\\n\\[\e[1;31m\\]\u2514\u2500\\[\e[1;33m\\]$ \\[\e[0m\\]'


# Modified Kali Linux prompt with arrow unicode
# Without hostname
#PS1=$'\[\e[1;31m\]\u250c\u2500\u2500\[\e[1;31m\](\[\e[1;33m\]\u\[\e[1;31m\])\u2500\[\e[1;31m\][\[\e[0;94m\]\w\[\e[1;31m\]]\n\[\e[1;31m\]\u2514\u2500\[\e[1;31m\]\u27a4 \[\e[0m\]'

# With hostname
PS1=$'\[\e[1;31m\]\u250c\u2500\u2500\[\e[1;31m\](\[\e[1;33m\]\u@\h\[\e[1;31m\])\u2500\[\e[1;31m\][\[\e[0;94m\]\w\[\e[1;31m\]]\n\[\e[1;31m\]\u2514\u2500\[\e[1;31m\]\u27a4 \[\e[0m\]'
