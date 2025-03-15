# Prints news headlines from RSS_FEED. 
#------------------------------------------------------
# INPUTS
#--------
# RSS_FEED: HTML link to the RSS Feed
# STR[1][2]: String headings printed above headlines
# TEST[1][2]: Text and brackground colours of STR[1][2]
#------------------------------------------------------
#!/bin/bash
#------------------------------------------------------
RSS_FEED=http://feeds.bbci.co.uk/news/world/rss.xml

STR1="BBC"
STR2="NEWS"
TEXT1='\033[47;2;31m'
TEXT2='\033[41;2;39m'
#------------------------------------------------------
NOCOLOR='\033[0m' 

# Print heading
echo -e "${TEXT1} ${STR1} ${TEXT2} ${STR2} ${NOCOLOR}"

# Grab and filter RSS feed
 curl --silent ${RSS_FEED} |  grep -E '(title>|description>)' | \
  head | sed -n '4,$p' | sed -e ' s/<title>//' -e 's/<description>//' \
  -e 's/<!\[CDATA\[//' | \sed -e 's/\]\]><\/title>//' | \
  sed -e 's/\]\]><\/description>//' 
