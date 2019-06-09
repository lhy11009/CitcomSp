import subprocess #module for running bashshell

#------------------------------------------------------
#lhy 20180530 dealing with images using the bash convert command
#n = 0 : convert filein to fileout, else append filein to fileout
#if has assigned size, then resize images
#-----------------------------------------------------
def convert_image(filein,fileout,n,oper='h',row=None,column=None):
    assert ((oper == 'h') or (oper == 'v')), "Error: operation for convert_image is either 'h' or 'v'" 
    if n == 0:
        bashCommand = 'convert %s' % (filein)    
    else:
    #either combine vertically or horizontally
        if oper == 'h':
            bashCommand = 'convert %s %s +append' % (fileout,filein)
        else:
            bashCommand = 'convert %s %s -append' % (fileout,filein)
    #resize picture
    if (row is not None) and (column is not None):
        appdix = '-resize %dx%d' % (row,column)
        bashCommand = ' '.join([bashCommand,appdix])
    bashCommand = ' '.join([bashCommand,fileout])
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return
