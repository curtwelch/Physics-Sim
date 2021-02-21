"""
    module crt
    ANSI CRT (VT100) Escape Codes

    functions all use print() to output codes to stdout.

    Only a very few basic function included.
    See: https://en.wikipedia.org/wiki/ANSI_escape_code
    for detailed documentation and history of these codes.

    2021-2-30 Curt Welch copied from other stuff I've used for years.
"""

from typing import Union


def clear():
    """ Clear screen without moving cursor. """
    print("\033[2J", end='')


def clear_and_home():
    """ Clear screen and move cursor to upper left corner. """
    clear()
    home()


def home():
    """
        Move cursor to top left corner.
        Alias for goto() or goto(1, 1)
    """
    goto(1, 1)


def goto(row: int = 1, col: int = 1):
    """ Position cursor to row, col (defaults to upper left corner)
        1,1 is top left corner.
        Rows move down lines.
        Columns move right across screen.
    """
    print(f"\033[{row};{col}H", end='')


def mode(m: Union[str, int] = 0):
    """
        Set video mode for characters.
        Defaults to normal mode.
        m == 0 or "normal" -- all effects off
        m == 1 or "bold"
        m == 2 or "underline"
        m == 5 or "blinking"
        m == 7 or "reverse video"
        For other values of m see:
        https://en.wikipedia.org/wiki/ANSI_escape_code#SGR
    """
    if m == "normal":
        m = 0
    elif m == "bold":
        m = 1
    elif m == "underline":
        m = 4
    elif m == "blinking":
        m = 5
    elif m == "reverse video":
        m = 7

    print(f"\033[{m}m", end='')


def unit_test():
    """ Simple unit testing of crt module code and of a terminal. """
    print()
    print("crt module unit test.")
    print()
    input("Press Enter to crt.clear() without moving cursor: ")
    clear()

    input("Press Enter to crt.clear_and_home(): ")
    clear_and_home()

    # print("01 Top left row of screen (row=0 column=0)")
    for i in range(1, 10):
        print(f"Line {i:02d}")

    print()
    input("Press Enter to crt.home() without clear: ")
    home()
    print("THIS IS HOME! ", end='')

    input("Press Enter for crt.goto() test.")
    clear()
    for row in range(1, 10):
        col = row*2 - 1
        goto(row, col)
        print(f"X crt.goto(row={row}, col={col})  ")
    print()

    input("Press Enter for crt.mode() test:")
    print()
    unit_test_mode("normal")
    unit_test_mode("bold")
    unit_test_mode("underline")
    unit_test_mode("blinking")
    unit_test_mode("reverse video")
    for m in range(10):
        unit_test_mode(m)

    print("Some terminals support mode values as high as 100")
    print()
    print("module crt unit test DONE!")


def unit_test_mode(m):
    print(f"This is ", end='')
    mode(m)
    print(f"crt.mode({m!r}) Testing 123...")
    mode()


if __name__ == '__main__':
    unit_test()
