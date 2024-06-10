class colors:
    ERROR = '\033[91m' #red
    PASS = '\033[92m' #green
    WARNING = '\033[93m' #yellow
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    END = '\033[0m'



def isValid(message):
    assert isinstance(message, str), f"{colors.ERROR} [ERROR] Lera's Trigger Contact Tool : Message {message} must be a string, but it is {type(message)} {colors.END}"
    return True 

def info(message):
   isValid(message)
   print("{}[INFO] Lera's Trigger Contact Tool : {}{}".format(colors.BLUE, colors.END, message))

def warning(message):
    isValid(message)
    print("{}[WARNING] Lera's Trigger Contact Tool : {}{}".format(colors.WARNING, colors.END, message))

def error(message):
    isValid(message)
    print("{}[ERROR] Lera's Trigger Contact Tool : {}{}".format(colors.ERROR, message, colors.END))


