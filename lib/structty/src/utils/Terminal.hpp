#pragma once
#include <unistd.h>

namespace Terminal {

struct Size {
    int rows;
    int cols;
};

struct MouseEvent {
    int  x;        // 0-based column
    int  y;        // 0-based row
    bool pressed;  // true = button down / motion, false = release
    bool moved;    // true = motion event (btn & 32 set in SGR)
};

// KEY_* 상수 — ncurses KEY_* 대응
constexpr int KEY_UP    = 0x100;
constexpr int KEY_DOWN  = 0x101;
constexpr int KEY_LEFT  = 0x102;
constexpr int KEY_RIGHT = 0x103;
constexpr int KEY_MOUSE = 0x104;

// raw mode 진입/복원 (SIGTERM/SIGINT atexit 핸들러 포함)
void enter_raw_mode();
void exit_raw_mode();

// 터미널 크기 (ioctl TIOCGWINSZ, 실패 시 80×24)
Size get_size();

// 화면 제어 (ANSI escape)
void clear();                       // \033[2J\033[H + fflush
void cursor_home();                 // \033[H — 커서만 (1,1)로, 화면 지우지 않음
void move_cursor(int row, int col); // \033[row+1;col+1H  (0-based 입력)
void clear_to_eol();                // \033[K

// 입력 버퍼 비우기
void flush_input();                 // tcflush(STDIN_FILENO, TCIFLUSH)

// 키 읽기 (blocking)
// 반환값: ASCII (1–127), KEY_* 상수, 또는 0 (unrecognized / error)
int  read_key();

// read_key()가 KEY_MOUSE를 반환한 직후 호출하여 마우스 좌표 취득
bool read_mouse(MouseEvent& ev);

} // namespace Terminal
