#!/bin/bash
#-------------------------------------------------------
# Converts user input Gregorian Calendar date (DDMMYYYY) 
# into Julian date (JD) and vice-versa
# ------------------------------------------------------
# INPUT
# USER_DATE_INPUT: Date from standard input. It must be
#                  in the format DD/MM/YYYY
# USER_TIME_INPUT: Time from standard input. It must be
#                  in the format HH:MM:SS. If time is
#                  not input (or input incorrectly), the
#                  JD is calculated at 00:00:00 for 
#                  USER_DATE_INPUT
#--------------------------------------------------------
# OUTPUT
# JD: Float of Julian Date on the Gregorian Calendar 
#     date specified by the user.
#========================================================
extract_var () {	
   var=`echo $1 | cut -d $2 -f$3`
   var=${var#0}
   echo "$var"
} 

num_check () {	
   if ! [[ $1 =~ ^[0-9]+$ ]] ; then
      echo "Error: non-numeric input" ;  exit 1
   fi
}
#--------------------------------------------------------
echo "Enter a Gregorian calendar date (dd/mm/yyyy):"
read user_date_input
#--------------------------------------------------------
# INITIAL DATE INPUT CHECK 
echo $user_date_input
echo $user_date_input | grep -q "/" 
if [[ $? != 0 ]] ; then 
 echo "Error: check date input format (dd/mm/yyyy)"
 exit 1
fi
#---------------------------------------------------------
d=$(extract_var "$user_date_input" '/' 1)
m=$(extract_var "$user_date_input" '/' 2)
yyyy=$(extract_var "$user_date_input" '/' 3)
#---------------------------------------------------------
# CHECK VALID DATE INPUT
# zero check
if [ -z $d ] || [ -z $m ] || [ -z $yyyy ] ; then
  echo "Error: invalid date"; exit 1
fi

num_check $d
num_check $m
num_check $yyyy

# Check date
if  [ $d -gt 31 ] || [ $m -gt 12 ] ; then
  echo "Error: invalid date"
  exit 1
fi

md=(4 6 9 11)
for (( i=0; i<${#md[@]}; i++));
 do
  if [ $m == ${md[$i]} ] && [ $d == 31 ] ; then
    echo "Error: max 30 days in your input month" ; exit 1
  fi
 done
# Leap year check for days in Feb
ly=$((yyyy%4))
if [ $ly != 0 ] && [ $m == 2 ] && [ $d == 29 ] ; then
  echo "Error: invalid day for a non-leap year" ; exit 1
fi
#-------------------------------------------------------- 
echo "Enter a time (hh:mm:ss):"
read user_time_input
#--------------------------------------------------------
# INITIAL TIME INPUT CHECK
echo $user_time_input | grep -q ":" 
if [[ $? != 0 ]] || [ -z $user_time_input ]; then 
 echo "Error: check time input format (hh:mm:ss)"
 echo "setting time to 00:00:00"
 user_time_input=00:00:00
fi
#-------------------------------------------------------
hr=$(extract_var "$user_time_input" ':' 1)
min=$(extract_var "$user_time_input" ':' 2)
sec=$(extract_var "$user_time_input" ':' 3)
#-------------------------------------------------------
# CHECK VALID TIME INPUT
num_check $hr
num_check $min
num_check $sec

if [ $hr -gt 23 ] || [ $min -gt 59 ] || [ $sec -gt 59 ] \
 ; then
  echo "Error: enter valid time"; exit 1
fi
#--------------------------------------------------------
# Correct inputs for Jan and Feb
if [ $m -le 2 ] ; then
  yyyy=$((yyyy-1))
  m=$m+12
fi
#--------------------------------------------------------
if [[ $yyyy -lt 100 ]] ; then
   i=0
else
   i=`echo "scale=3;$yyyy / 100" | bc | cut -d '.' -f1`
fi

if [[ $yyyy -lt 400 ]] ; then
   j=0
else
   j=`echo "scale=3;$i / 4" | bc | cut -d '.' -f1`
fi

k=`echo "scale=3;2 - $i + $j" | bc`
x=`echo 365.25\*"($yyyy+4716)" | bc | cut -d '.' -f1`
z=`echo 30.6001\*"($m+1)" | bc | cut -d '.' -f1`

sum=`echo "$k+$d+$x+$z" | bc | cut -d '.' -f1`  

JD=$sum-1524.5
#-------------------------------------------------------
# Time conversion
min_hr=`echo "scale=4; $min / (60 * 24)" | bc`
sec_hr=`echo "scale=5; $sec / (60 * 60 * 24)" | bc`
time=`echo "scale=5;$hr / 24 + $min_hr + $sec_hr" | bc`
string_output="Julian date is "
JD=`echo "$JD +$time" | bc`
string_output+="$JD"
echo $string_output
