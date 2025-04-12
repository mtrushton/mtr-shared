 wget randomfunfacts.com -qO- | sed -n "/<strong>/{s;^.*<i>\(.*\)</i>.*$;\1;p}"
