#include "Terminal.hpp"
#include <termios.h>
#include <sys/ioctl.h>
#include <sys/select.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <cerrno>

namespace Terminal {

// --- internal state ---
static struct termios g_saved_termios;
static bool           g_raw_active = false;

// Buffered mouse event from the most recent SGR sequence
static int  g_mouse_x       = 0;
static int  g_mouse_y       = 0;
static bool g_mouse_pressed = false;
static bool g_mouse_moved   = false;

// --- cleanup (async-signal-safe: only write() + tcsetattr()) ---
static void do_exit_raw() {
    if (!g_raw_active) return;
    g_raw_active = false;

    // Exit alternate screen, disable mouse tracking, restore cursor in one write
    static const char disable[] =
        "\033[?1049l"   // exit alternate screen buffer (restores original content)
        "\033[?1003l"   // any-event off
        "\033[?1002l"   // button-event off
        "\033[?1000l"   // mouse off
        "\033[?1006l"   // SGR format off
        "\033[?25h";    // show cursor
    write(STDOUT_FILENO, disable, sizeof(disable) - 1);

    tcsetattr(STDIN_FILENO, TCSAFLUSH, &g_saved_termios);
}

static void atexit_handler()  { do_exit_raw(); }
static void sig_handler(int)  { do_exit_raw(); _exit(0); }

// --- public API ---

void enter_raw_mode() {
    if (g_raw_active) return;

    if (tcgetattr(STDIN_FILENO, &g_saved_termios) == -1) return;
    g_raw_active = true;

    struct termios raw = g_saved_termios;
    // Disable canonical mode, echo, signals, flow control
    raw.c_iflag &= ~(unsigned)(BRKINT | ICRNL | INPCK | ISTRIP | IXON);
    raw.c_oflag &= ~(unsigned)(OPOST);
    raw.c_cflag |=  (unsigned)(CS8);
    raw.c_lflag &= ~(unsigned)(ECHO | ICANON | IEXTEN | ISIG);
    raw.c_cc[VMIN]  = 1;  // block until at least 1 byte
    raw.c_cc[VTIME] = 0;  // no timeout
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);

    // Enable mouse tracking, enter alternate screen, hide cursor
    static const char enable[] =
        "\033[?1000h"   // mouse reporting on
        "\033[?1002h"   // button-event tracking
        "\033[?1003h"   // any-event (all mouse movements)
        "\033[?1006h"   // SGR extended format (supports large coordinates)
        "\033[?1049h"   // enter alternate screen buffer (clean slate, restores on exit)
        "\033[?25l";    // hide cursor
    write(STDOUT_FILENO, enable, sizeof(enable) - 1);
    fflush(stdout);

    atexit(atexit_handler);
    signal(SIGTERM, sig_handler);
    signal(SIGINT,  sig_handler);
    // SIGWINCH default (Ignore) is sufficient; get_size() is called each frame
}

void exit_raw_mode() {
    do_exit_raw();
}

Size get_size() {
    struct winsize ws;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == -1 || ws.ws_row == 0) {
        return {24, 80};  // safe fallback
    }
    return {(int)ws.ws_row, (int)ws.ws_col};
}

void clear() {
    // Erase display then move cursor to home
    write(STDOUT_FILENO, "\033[2J\033[H", 7);
    fflush(stdout);
}

void cursor_home() {
    // Move cursor to (1,1) without erasing — used for flicker-free per-frame redraw
    write(STDOUT_FILENO, "\033[H", 3);
}

void move_cursor(int row, int col) {
    char buf[24];
    int n = snprintf(buf, sizeof(buf), "\033[%d;%dH", row + 1, col + 1);
    if (n > 0) fwrite(buf, 1, (size_t)n, stdout);
}

void clear_to_eol() {
    fwrite("\033[K", 1, 3, stdout);
}

void flush_input() {
    tcflush(STDIN_FILENO, TCIFLUSH);
}

int read_key() {
    unsigned char c = 0;
    ssize_t n = read(STDIN_FILENO, &c, 1);
    if (n <= 0) return 0;

    if (c != 0x1B) return (int)c;  // plain ASCII

    // ESC received — check for sequence continuation within 50 ms
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);
    struct timeval tv = {0, 50000};
    if (select(STDIN_FILENO + 1, &fds, nullptr, nullptr, &tv) <= 0) {
        return 0x1B;  // bare ESC (not bound to any action, but safe to return)
    }

    // Expect '['  (CSI)
    unsigned char c2 = 0;
    if (read(STDIN_FILENO, &c2, 1) <= 0) return 0x1B;
    if (c2 != '[') return 0;

    // Read next byte to dispatch
    unsigned char c3 = 0;
    if (read(STDIN_FILENO, &c3, 1) <= 0) return 0;

    if (c3 == '<') {
        // SGR mouse sequence: \033[<btn;x;yM or \033[<btn;x;ym
        // Read digits and delimiters until 'M' or 'm'
        char numstr[32];
        int  nlen = 0;
        while (nlen < (int)(sizeof(numstr) - 1)) {
            unsigned char ch = 0;
            if (read(STDIN_FILENO, &ch, 1) <= 0) break;
            if (ch == 'M' || ch == 'm') {
                numstr[nlen] = '\0';
                int btn = 0, x = 0, y = 0;
                sscanf(numstr, "%d;%d;%d", &btn, &x, &y);
                g_mouse_x       = x - 1;          // 1-based → 0-based
                g_mouse_y       = y - 1;
                g_mouse_pressed = (ch == 'M');
                g_mouse_moved   = (btn & 32) != 0; // bit 5 = motion
                return KEY_MOUSE;
            }
            numstr[nlen++] = (char)ch;
        }
        return 0;
    }

    // Cursor-key sequences: \033[A/B/C/D
    switch (c3) {
        case 'A': return KEY_UP;
        case 'B': return KEY_DOWN;
        case 'C': return KEY_RIGHT;
        case 'D': return KEY_LEFT;
    }

    // Discard remaining bytes of unrecognised multi-char sequence
    // (e.g. \033[3~ for Delete) — ignore, drain with a non-blocking read
    {
        struct timeval drain_tv = {0, 0};
        fd_set dfds;
        FD_ZERO(&dfds);
        FD_SET(STDIN_FILENO, &dfds);
        while (select(STDIN_FILENO + 1, &dfds, nullptr, nullptr, &drain_tv) > 0) {
            unsigned char dummy;
            if (read(STDIN_FILENO, &dummy, 1) <= 0) break;
            if (dummy >= 0x40 && dummy <= 0x7E) break;  // final byte reached
            FD_ZERO(&dfds);
            FD_SET(STDIN_FILENO, &dfds);
            drain_tv = {0, 0};
        }
    }

    return 0;  // unrecognised — ignore
}

bool read_mouse(MouseEvent& ev) {
    ev.x       = g_mouse_x;
    ev.y       = g_mouse_y;
    ev.pressed = g_mouse_pressed;
    ev.moved   = g_mouse_moved;
    return true;
}

} // namespace Terminal
