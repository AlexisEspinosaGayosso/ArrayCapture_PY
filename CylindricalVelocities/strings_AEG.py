#AEG: This module is intended to keep all the functions used for manipulating or processing strings



################################################################
#is_number(s)

#----------
#2017.04.11
#Function to evaluate if a string is a number.

#Copied-Pasted from website: http://pythoncentral.io/how-to-check-if-a-string-is-a-number-in-python-including-unicode/
#It also works with non-ASCII numbers encoded in unicode
#----------

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
 
    finally:
        return False
################################################################
