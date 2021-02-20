"""
    module crt
    ANSI CRT (VT100) Control

    2021-2-30 Curt Welch copied from other stuff I've used for years.
"""


def clear():
    """ Clear screen without moving cursor, """
    print("\033[2J", end='')


def clear_and_home(row: int = 0, col: int = 0):
    """ Clear and move cursor (upper left corner by default) """
    clear()
    goto(row, col)


def goto(row: int = 0, col: int = 0):
    """ Move cursor to row, col (defaults to upper left corner)
        0,0 is top left corner.
        Rows move down lines.
        Columns move right across characters.
    """
    print(f"\033[{row};{col}H", end='')


def mode(mode_str: str = None, m_value: int = None):
    # TODO finish writing and testing this
    """
        Set video mode for characters drawn.
        Defaults to normal mode -- all special effects turned off.
        mode_str == NONE -> normal mode -- all off
        mode_str == "normal" same as NONE
        mode_str == "bold"
        mode_str == "underline"
        mode_str == "blinking"
        mode_str == "reverse video"
    """
    if m_value is not None:
        if mode_str is None:
            m_value = 0  # normal mode
        elif mode_str == "normal":
            m_value = 0
        elif mode_str == "bold":
            m_value = 1
        elif mode_str == "underline":
            m_value = 2
        elif mode_str == "blinking":
            m_value = 5
        elif mode_str == "reverse video":
            m_value = 7

    if m_value not in [0, 1, 2, 5, 7]:
        raise ValueError("Invalid mode")

    print(f"\033[{m_value}m", end='')
