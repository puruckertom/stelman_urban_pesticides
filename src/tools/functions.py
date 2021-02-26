import logging, os, sys, re, __main__ as main
from datetime import datetime
# from paths import * # this causes errors, so just paste the part we need
main_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# set up logger
## borrowed from source: "https://docs.python.org/2/howto/logging.html#logging-basic-tutorial"
logging.basicConfig(filename='pyabc.log', level=logging.INFO)

'''
Sets up the function 'loginfo' with the appropriate prefix for the current script
 Input: script <str> -Script's ID, will be used in prefix of logging messages (ex. '01a', '08', '10')-
 Output: loginfo <function> -Function that takes in a message and logs it with the appropriate prefix for current script-
'''
def log_prefixer(script):
    def loginfo(text):
        logging.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ":: " + script + ": " + text)
    def logerror(text):
        logging.error(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ":: " + script + ": " + text)
    return(loginfo, logerror)

# for 01a and 03
# Borrowed from "https://www.oreilly.com/library/view/python-cookbook/0596001673/ch04s16.html"
'''
Breaks a directory of file path into all its parts
 Input: path <str> -Path (must be valid) to be split up-
 Output: <list of str> -a list of the path's components in ancestral order-
'''
def splitall(path):
    allparts = []
    while 1:
        parts = os.path.split(path)
        if parts[0] == path:  # sentinel for absolute paths
            if path not in ["/","\\","\\\\"]:
                allparts.insert(0, parts[0])
            else:
                allparts.insert(0, "")
            break
        elif parts[1] == path: # sentinel for relative paths
            if path not in ["/","\\","\\\\"]:
                allparts.insert(0, parts[1])
            else:
                allparts.insert(0, "")
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])
    return allparts

# for 01a and 03
'''
Corrects absolute paths in input file, so the input files work on any computer, not just the author's
 Inputs: inp_path <str> -Path to input file (Optional if filelines provided)-
   filelines <list of str> -Lines of the file to clean up (Optional if inp_path provided. Dominates inp_path)-
   new_path <str> -Path to file to be overwritten with corrected content (Optional*)-
   Log <Bool> -Indicator for whether to write logging messages to log file (True/default), or to print them instead (False)-
 Output: *ONLY IF new_path missing: <list of str> The cleaned up lines of the file
'''
def replace_infile_abspaths(inp_path = None, filelines = None, new_path = None, Log = True):
    # raise TypeError if both inp_path and filelines have not been provided
    if not (inp_path or filelines):
        raise TypeError("Missing arguments for both 'inp_path' and 'filelines' (1 of the 2 must be provided)")
    # if filelines is not provided, extract them from the inp_file manually
    if not filelines:
        # Set up logging for inner-file operations
        try:
            script = script
        except NameError:
            try:
                script = main.script
            except AttributeError:
                script = os.path.basename(main.__file__)[:3]
        # Write info to log file if Log arg is True, or print info if False
        if Log:
            try:
                loginfo = loginfo
            except NameError:
                loginfo, logerror = log_prefixer(script)
        else:
            loginfo = print

        # read the input file and extract its contents as a list
        loginfo("Opening file <" + inp_path + "> to read content out of.")
        ip_file = open(inp_path, 'r')
        filelines = ip_file.readlines()
        loginfo("Closing file <" + inp_path + ">.")
        ip_file.close()
    
    # the first absolute path to correct, listified
    path1cols = filelines[50].split()
    # remember, there might be a space in the filepath, meaning that the split function could have created two elements, not 1
    # so instead, make a new list using the first five, a space holder, and the last two elements of the original list
    path1cols = path1cols[:5] + [""] + path1cols[-2:]
    # the corrected element of the listified line
    path1cols[5] = '"' + os.path.join(main_path, "master", "weather", "swmm_wet.txt") + '"'
    # insert the correction and unlistify!
    filelines[50] = "\t".join(path1cols) + "\n"

    # the second absolute path to correct, listified
    path2cols = filelines[1384].split()
    # remember, there might be a space in the filepath, meaning that the split function could have created two elements, not 1
    # so instead, make a new list using the first 2 elements of the original list and a space holder
    path2cols = path2cols[:2] + [""]
    # the corrected element of the listified line
    path2cols[2] = '"' + os.path.join(main_path, "master", "app_rate_output_for_swmm_48rain.txt") + '"'
    # insert the correction and unlistify!
    filelines[1384] = "\t".join(path2cols) + "\n"

    # the third absolute path to correct, listified
    path3cols = filelines[9306].split()
    # remember, there might be a space in the filepath, meaning that the split function could have created two elements, not 1
    # so instead, make a new list using the first element of the original list and a space holder
    path3cols = path3cols[:1] + [""]
    # the corrected element of the listified line
    path3cols[1] = '"' + os.path.join(main_path, "master", "nplesant.jpg") + '"'
    # insert the correction and unlistify!
    filelines[9306] = "\t".join(path3cols) + "\n"

    if new_path:
        # copy, write out file
        loginfo("Opening file <" + new_path + "> to overwrite with edited content.")
        new_file = open(new_path, 'w')
        new_file.writelines(filelines)
        loginfo("Closing file <" + new_path + ">.")
        new_file.close()
    else: 
        return filelines

# for 01b and 05
'''
Saves data frame to specified .csv file and returns it
 Inputs: df <pandas.DataFrame> -Data frame to export to .csv-
   csv <str> -Path to .csv file where data frame is to be exported-
   msg <Bool or str> 
     -If str, message body to write to logging file-
     -If bool, indicator of whether to write default message to logging file (True/default), or not to log at all (False)-
 Output: df <pandas.DataFrame> -Same data frame from input-
'''
def save_and_continue(df,csv,msg = True):
    if not isinstance (msg,str):
        if msg == True:
            bn = os.path.basename(csv)
            dn = os.path.basename(os.path.dirname(csv))
            msg = "Saving intermediate version of data to <" + bn + "> in <" + dn + ">."
    if msg:
        # Set up logging for inner-file operations
        try: 
            loginfo = loginfo
        except NameError:
            try:
                script = script
            except NameError:
                try:
                    script = main.script
                except AttributeError:
                    script = os.path.basename(main.__file__)[:3]
            loginfo, logerror = log_prefixer(script) 
        loginfo(msg)
    df.to_csv(csv)
    return(df)

# for 01b, 04, and 05
'''
Saves data frame to specified .csv file and returns message of completion
 Inputs: df <pandas.DataFrame> -Data frame to export to .csv-
   csv <str> -Path to .csv file where data frame is to be exported-
   msg <Bool or str> 
     -If str, message body to write to logging file-
     -If bool, indicator of whether to write default message to logging file (True/default), or not to log at all (False)-
 Output: <str> -Message of completion-
'''
def save_and_finish(df,csv,msg = True):
    if not isinstance (msg,str):
        if msg == True:
            bn = os.path.basename(csv)
            dn = os.path.basename(os.path.dirname(csv))
            msg = "Saving final version of data to <" + bn + "> in <" + dn + ">."
    if msg:
        # Set up logging for inner-file operations
        try: 
            loginfo = loginfo
        except NameError:
            try:
                script = script
            except NameError:
                try:
                    script = main.script
                except AttributeError:
                    script = os.path.basename(main.__file__)[:3]
            loginfo, logerror = log_prefixer(script) 
        loginfo(msg)
    df.to_csv(csv)
    return("Finished " + dn)


# for 03

'''
Helper function of edit_file:
Create new line of .inp file with lhs simulated versions of values and text in old .inp file
 Inputs: fileline <str> -line from original .inp file-
   Col <int> -column in the .inp file that contains the info you want to preserve and clean- 
   sim <float> -the simulated value to repace the old (observed) one with-
 Output: newline <str> -custom version of line for new file-
'''
def edit1line(fileline, Col, sim):
    listline = fileline.split()
    listline[Col] = str(sim)
    newline = ' '.join([str(item) for item in listline]) + "\n"
    return(newline)

# for 03

'''
Edit the new file to be cleaner
 Inputs: Ite <int> -index of current simulation-
   Num <int> -Number of rows to clean-
   row_0 <int> -Number of rows to skip-
   parameter <str> -name of lhs design parameter (will become name of column in the .csv  file)-
   Col <int> -column in the .inp file that contains the info you want to preserve and clean- 
   flines <list of str> -lines of the file to clean up-
 Output: cleaned up lines of file given
'''
def editted_lines(swmm_dict, Num, row_0, parameter, Col, flines):
    # value of "parameter" in current lhs simulation
    sim = swmm_dict[parameter]
    # print(sim)
    return([edit1line(flines[row_0 + i], Col, sim) for i in range(Num)])

'''
Delete the file, files, or directory at a given path location
 Inputs: path <path (str)> -path to the item to be deleted-
   *args <path (str)> -path to additional items to be deleted-
 (No output)
'''
def rm(path, *args):
    try:
        assert os.system('rm -r ' + path) == 0, "original path removal attempt failed"
        print("original path removal successful")
    except AssertionError as err1:
        #print(err1)
        path = re.sub(r'\\', r'/', path)
        try:
            assert os.system('rm -r ' + path) == 0, "second path removal attempt failed"
            print("second path removal successful")
        except AssertionError as err2:
            print(err2)
    for p in args:
        rm(path = p)
