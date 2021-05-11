import os, sys, subprocess

runAnalyzer = True

signalFile = ''
signalName = 'Signal'

backgroundFiles = ['root://eoscms.cern.ch://eos/cms/store/user/lzygala/HVV/ttHTobb/4A5BE1BE-DA13-E811-BEB0-AC1F6B1AEFFC.root']
backgroundNames = ['']

def argConversion(arg):
    if type(arg) in [float, int]:
        return str(arg)
    elif type(arg) in [str]:
        return '"' + arg + '"'
    else:
        return '"' + str(arg) + '"'

def argConnection(arglist=[]):
    if not (arglist is [] or arglist is None):
        shellargv = list(arglist) 
        arglist_conv = list()
        for arg in arglist:
            arglist_conv.append(argConversion(arg))
        return '('+','.join(arglist_conv)+')'
    else:
        return ''

def runMacro(macroName, arglist=None, splash=False, interprete=False, batch=True):
    shellCommand = ['root']
    if interprete is False:
        shellCommand.append("-q")
    if splash is False:
        shellCommand.append("-l")
    if batch is True:
        shellCommand.append("-b")
    shellCommand.append(macroName+argConnection(arglist))
    print("Run Macro", shellCommand)
    a = subprocess.Popen(shellCommand)
    return a

def main():
    analyzerArgList = [signalFile,
                       signalName,
                       backgroundFiles,
                       backgroundNames]

    if runAnalyzer:
        print("Running Analyzer")
        analyzer = runMacro('Analyzers/EnvelopeAnalyzer.C++', analyzerArgList)
        analyzer.wait()


if __name__ == "__main__":
    main()