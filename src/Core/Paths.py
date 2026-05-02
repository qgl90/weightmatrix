import os 
def expand_env_var(path: str) -> str:
    """
    Expand environment variables in a path.
    Supports both Unix-style ($VAR) and Windows-style (%VAR%) variables.
    
    Examples:
        "~/docs"          → "/home/user/docs"
        "$HOME/docs"      → "/home/user/docs"
        "%USERPROFILE%\\data" → "C:\\Users\\username\\data"
    """
    if not path:
        return path
    
    # os.path.expandvars handles both $VAR and %VAR% formats
    expanded = os.path.expandvars(path)
    
    # Also expand ~ (home directory)
    expanded = os.path.expanduser(expanded)
    
    return expanded