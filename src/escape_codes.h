/* **************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke   *
 *   mihrke@uni-goettingen.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/* Curtesy to Oliver Bendix
** tel 0551-5176-431, fax 0551-5176-409, mail oliver[at]nld.ds.mpg.de **
**  Nonlinear Dynamics Group, MPI for Dynamics And Self-Organization  **
**         Bunsenstrasse 10, 37073 Goettingen, Postfach 2853          **
**                   http://www.nld.ds.mpg.de                         **
*/

/** \file
* \brief \ref status_stable Definitions of print features.

* If the shell environment supports this, there are special features
* we can use. We also included ANSI CSI (Control Sequences Intitiator)
* codes.
* @note This works only, if WITH_ANSI_ESCAPE_CODES was defined before
* (see in Makefile.$PLATFORM)
*/
#ifndef ESCAPE_CODES_H
#define ESCAPE_CODES_H

#define ESCAPE_RESET        0
#define ESCAPE_BRIGHT         1
#define ESCAPE_BOLD         1
#define ESCAPE_DIM        2
#define ESCAPE_FAINT        2
#define ESCAPE_ITALIC        3
#define ESCAPE_UNDERLINE     4
#define ESCAPE_BLINK        5
#define ESCAPE_BLINKRAPID    6
#define ESCAPE_REVERSE        7
#define ESCAPE_HIDDEN        8
#define ESCAPE_CONCEAL        8
#define ESCAPE_DUNDERLINE     21

#define ESCAPE_BRIGHT_OFF        22
#define ESCAPE_BOLD_OFF          22
#define ESCAPE_DIM_OFF           22
#define ESCAPE_FAINT_OFF         22
#define ESCAPE_ITALIC_OFF        23
#define ESCAPE_UNDERLINE_OFF     24
#define ESCAPE_BLINK_OFF         25
#define ESCAPE_REVERSE_OFF       27
#define ESCAPE_HIDDEN_OFF        28
#define ESCAPE_CONCEAL_OFF       29

#define ESCAPE_BLACK         30
#define ESCAPE_RED        31
#define ESCAPE_GREEN        32
#define ESCAPE_YELLOW        33
#define ESCAPE_BLUE        34
#define ESCAPE_MAGENTA        35
#define ESCAPE_CYAN        36
#define ESCAPE_WHITE        37

#define ESCAPE_GREY         90
#define ESCAPE_LRED        91
#define ESCAPE_LGREEN        92
#define ESCAPE_LYELLOW        93
#define ESCAPE_LBLUE        94
#define ESCAPE_LMAGENTA        95
#define ESCAPE_LCYAN        96
#define ESCAPE_LWHITE        97

#if defined (_WITH_ANSI_ESCAPE_CODES)

 #define ESCAPE_RESET_STR       "\033[0m"
 #define ESCAPE_BRIGHT_STR       "\033[1m"
 #define ESCAPE_BOLD_STR           "\033[1m"
 #define ESCAPE_DIM_STR           "\033[2m"
 #define ESCAPE_FAINT_STR          "\033[2m"
 #define ESCAPE_UNDERLINE_STR        "\033[4m"
 #define ESCAPE_BLINK_STR       "\033[5m"
 #define ESCAPE_BLINKRAPID_STR       "\033[6m"
 #define ESCAPE_REVERSE_STR       "\033[7m"
 #define ESCAPE_HIDDEN_STR       "\033[8m"
 #define ESCAPE_CONCEAL_STR       "\033[8m"
 #define ESCAPE_DUNDERLINE_STR        "\033[21m"

 #define ESCAPE_BRIGHT_OFF_STR     "\033[22m"
 #define ESCAPE_BOLD_OFF_STR       "\033[22m"
 #define ESCAPE_DIM_OFF_STR        "\033[22m"
 #define ESCAPE_FAINT_OFF_STR      "\033[22m"
 #define ESCAPE_ITALIC_OFF_STR     "\033[23m"
 #define ESCAPE_UNDERLINE_OFF_STR  "\033[24m"
 #define ESCAPE_BLINK_OFF_STR      "\033[25m"
 #define ESCAPE_REVERSE_OFF_STR    "\033[27m"
 #define ESCAPE_HIDDEN_OFF_STR     "\033[28m"
 #define ESCAPE_CONCEAL_OFF_STR    "\033[29m"

 #define ESCAPE_FGBLACK_STR       "\033[30m"
 #define ESCAPE_FGRED_STR       "\033[31m"
 #define ESCAPE_FGGREEN_STR       "\033[32m"
 #define ESCAPE_FGYELLOW_STR       "\033[33m"
 #define ESCAPE_FGBLUE_STR       "\033[34m"
 #define ESCAPE_FGMAGENTA_STR       "\033[35m"
 #define ESCAPE_FGCYAN_STR       "\033[36m"
 #define ESCAPE_FGWHITE_STR       "\033[37m"
 #define ESCAPE_FGGREY_STR       "\033[90m"
 #define ESCAPE_FGLRED_STR       "\033[91m"
 #define ESCAPE_FGLGREEN_STR       "\033[92m"
 #define ESCAPE_FGLYELLOW_STR       "\033[93m"
 #define ESCAPE_FGLBLUE_STR       "\033[94m"
 #define ESCAPE_FGLMAGENTA_STR       "\033[95m"
 #define ESCAPE_FGLCYAN_STR       "\033[96m"
 #define ESCAPE_FGLWHITE_STR       "\033[97m"

 #define ESCAPE_BGBLACK_STR       "\033[40m"
 #define ESCAPE_BGRED_STR       "\033[41m"
 #define ESCAPE_BGGREEN_STR       "\033[42m"
 #define ESCAPE_BGYELLOW_STR       "\033[43m"
 #define ESCAPE_BGBLUE_STR       "\033[44m"
 #define ESCAPE_BGMAGENTA_STR       "\033[45m"
 #define ESCAPE_BGCYAN_STR       "\033[46m"
 #define ESCAPE_BGWHITE_STR       "\033[47m"
 #define ESCAPE_BGGREY_STR       "\033[100m"
 #define ESCAPE_BGLRED_STR       "\033[101m"
 #define ESCAPE_BGLGREEN_STR       "\033[102m"
 #define ESCAPE_BGLYELLOW_STR       "\033[103m"
 #define ESCAPE_BGLBLUE_STR       "\033[104m"
 #define ESCAPE_BGLMAGENTA_STR       "\033[105m"
 #define ESCAPE_BGLCYAN_STR       "\033[106m"
 #define ESCAPE_BGLWHITE_STR       "\033[107m"

 /** Prints to file with ANSI attribute attr, foregroundcolor fg and backgroundcolor bg. */
 #define FPRINTFTEXTALL(file,attr,fg,bg)  fprintf(file,"%c[%d;%d;%dm",0x1b,attr,fg,bg+10)
 /** * @see FPRINTFTEXTALL(file,attr,fg,bg) , but prints to stdout */
 #define PRINTFTEXTALL(attr,fg,bg)        FPRINTFTEXTALL(stdout,attr,fg,bg)

 /** Prints to file with ANSI attribute attr */
 #define FPRINTFTEXTATTR(file,attr)       fprintf(file,"%c[%dm",0x1b,attr)
 /** @see FPRINTFTEXTATTR(file,attr), but prints to stdout */
 #define PRINTFTEXTATTR(attr)             FPRINTFTEXTATTR(stdout,attr)

 /** Prints to file with ANSI foregroundcolor fg. */
 #define FPRINTFTEXTFG(file,fg)           fprintf(file,"%c[%dm",fg)
 /** @see FPRINTFTEXTFG(file,fg), but prints to stdout */
 #define PRINTFTEXTFG(fg)                 FPRINTFTEXTFG(stdout,fg)

 /** Prints to file with ANSI backgroundcolor bg. */
 #define FPRINTFTEXTBG(file,bg)           fprintf(file,"%c[%dm",0x1b,bg+10)
 /** @see FPRINTFTEXTBG(file,bg), but prints to stdout */
 #define PRINTFTEXTBG(bg)                 FPRINTFTEXTALL(stdout,bg)

 /** reset all ANSI attributes in file */
 #define FPRINTFTEXTRESET(file)           fprintf(file,"%c[%dm",0x1b,ESCAPE_RESET)
 /** @see FPRINTFTEXTRESET(file), but print out to sdtout */
 #define PRINTFTEXTRESET()                FPRINTFTEXTRESET(stdout)

 /** Insert n Characters. Make room for n characters at current position */
 #define FPRINTFTEXT_ICH(file,n)       fprintf(file,"%c[%d@",0x1b,(n))
 #define PRINTFTEXT_ICH(n)                FPRINTFTEXT_ICH(stdout,n)

 /** Moves the cursor n cells in the given direction.
     If the cursor is already at the edge of the screen, this has no effect.
     CUU: Up; CUD: Down; CUF: Forward; CUB: Back */
 #define FPRINTFTEXT_CUU(file,n)       fprintf(file,"%c[%dA",0x1b,(n))
 #define FPRINTFTEXT_CUD(file,n)       fprintf(file,"%c[%dB",0x1b,(n))
 #define FPRINTFTEXT_CUF(file,n)       fprintf(file,"%c[%dC",0x1b,(n))
 #define FPRINTFTEXT_CUB(file,n)       fprintf(file,"%c[%dD",0x1b,(n))
 #define PRINTFTEXT_CUU(n)           FPRINTFTEXT_CUU(stdout,n)
 #define PRINTFTEXT_CUD(n)               FPRINTFTEXT_CUD(stdout,n)
 #define PRINTFTEXT_CUF(n)           FPRINTFTEXT_CUU(stdout,n)
 #define PRINTFTEXT_CUB(n)           FPRINTFTEXT_CUB(stdout,n)

 /** Moves cursor to beginning of the line n lines down. */
 #define FPRINTFTEXT_CNL(file,n)       fprintf(file,"%c[%dE",0x1b,(n))
 #define PRINTFTEXT_CNL(n)           FPRINTFTEXT_CNL(stdout,n)

 /** Moves cursor to beginning of the line n lines up. */
 #define FPRINTFTEXT_CPL(file,n)       fprintf(file,"%c[%dF",0x1b,(n))
 #define PRINTFTEXT_CPL(n)           FPRINTFTEXT_CPL(stdout,n)

 /** Moves the cursor to column n. */
 #define FPRINTFTEXT_CHA(file,n)       fprintf(file,"%c[%dG",0x1b,(n))
 #define PRINTFTEXT_CHA(n)           FPRINTFTEXT_CHA(stdout,n)

 /** Moves the cursor to row n, column m. */
 #define FPRINTFTEXT_CUP(file,n,m)       fprintf(file,"%c[%d;%dH",0x1b,(n),(m))
 #define PRINTFTEXT_CUP(n,m)           FPRINTFTEXT_CUP(stdout,n,m)

 /** Cursor Horizontal Tabulation for n tabs. */
 #define FPRINTFTEXT_CHT(file,n)       fprintf(file,"%c[%dI",0x1b,(n))
 #define PRINTFTEXT_CHT(n)           FPRINTFTEXT_CHT(stdout,n)

 /** Clears part of the screen. If n is zero, clear from cursor to end of
    screen. If n is one, clear from cursor to beginning of the screen.
    If n is two, clear entire screen. */
 #define FPRINTFTEXT_ED(file,n)       fprintf(file,"%c[%dJ",0x1b,(n))
 #define PRINTFTEXT_ED(n)           FPRINTFTEXT_ED(stdout,n)

 /** Erases part of the line. If n is zero, clear from cursor to the end of
    the line. If n is one, clear from cursor to beginning of the line. If
    n is two, clear entire line. Cursor position does not change. */
 #define FPRINTFTEXT_EL(file,n)       fprintf(file,"%c[%dK",0x1b,(n))
 #define PRINTFTEXT_EL(n)           FPRINTFTEXT_EL(stdout,n)

 /** Insert n lines before current line, scroll current line down. */
 #define FPRINTFTEXT_IL(file,n)       fprintf(file,"%c[%dL",0x1b,(n))
 #define PRINTFTEXT_IL(n)           FPRINTFTEXT_IL(stdout,n)

 /** Delete n lines from current line downward, do not scroll. */
 #define FPRINTFTEXT_DL(file,n)       fprintf(file,"%c[%dM",0x1b,(n))
 #define PRINTFTEXT_DL(n)           FPRINTFTEXT_DL(stdout,n)

 /** Delete Character n characters from the current cursor position. */
 #define FPRINTFTEXT_DCH(file,n)       fprintf(file,"%c[%dP",0x1b,(n))
 #define PRINTFTEXT_DCH(n)           FPRINTFTEXT_DCH(stdout,n)

 /** Scroll whole page up by n lines. New lines are added at the bottom. */
 #define FPRINTFTEXT_SU(file,n)       fprintf(file,"%c[%dS",0x1b,(a))
 #define PRINTFTEXT_SU(n)           FPRINTFTEXT_SU(stdout,n)

 /** Scroll whole page down by n lines. New lines are added at the top. */
 #define FPRINTFTEXT_SD(file,n)       fprintf(file,"%c[%dT",0x1b,(n))
 #define PRINTFTEXT_SD(n)           FPRINTFTEXT_SD(stdout,n)

 /** Moves the cursor to row n, column m. Same as CUP */
 #define FPRINTFTEXT_HVP(file,n,m)       fprintf(file,"%c[%d;%df",0x1b,(n),(m))
 #define PRINTFTEXT_HVP(n,m)           FPRINTFTEXT_HVP(stdout,n,m)

 /** Sets SGR (Select Graphic Rendition) parameters. */
 #define FPRINTFTEXT_SGR(file,n)       fprintf(file,"%c[%dm",0x1b,(n))
 #define PRINTFTEXT_SGR(n)           FPRINTFTEXT_SGR(stdout,n)

 /** Reports the cursor position to the application as (as though typed at
     the keyboard) ESC[n;mR, where n is the row and m is the column. */
 #define FPRINTFTEXT_DSR(file)              fprintf(file,"%c[6n",0x1b)
 #define PRINTFTEXT_DSR()           FPRINTFTEXT_DSR(stdout)

 /** Saves the cursor position. */
 #define FPRINTFTEXT_SCP(file)              fprintf(file,"%c[s",0x1b)
 #define PRINTFTEXT_SCP()           FPRINTFTEXT_SCP(stdout)

 /** Restores the cursor position. */
 #define FPRINTFTEXT_RCP(file)              fprintf(file,"%c[u",0x1b)
 #define PRINTFTEXT_RCP()           FPRINTFTEXT_RCP(stdout)

 /** Saves the cursor position and attributes. */
 #define FPRINTFTEXT_SCA(file)              fprintf(file,"%c7",0x1b)
 #define PRINTFTEXT_SCA()           FPRINTFTEXT_SCA(stdout)

 /** Restores the cursor position and attributes. */
 #define FPRINTFTEXT_RCA(file)              fprintf(file,"%c8",0x1b)
 #define PRINTFTEXT_RCA()           FPRINTFTEXT_RCA(stdout)

#else

 #define ESCAPE_RESET_STR          ""
 #define ESCAPE_BRIGHT_STR         ""
 #define ESCAPE_BOLD_STR           ""
 #define ESCAPE_DIM_STR            ""
 #define ESCAPE_FAINT_STR          ""
 #define ESCAPE_UNDERLINE_STR      ""
 #define ESCAPE_DUNDERLINE_STR     ""
 #define ESCAPE_BLINK_STR          ""
 #define ESCAPE_BLINKRAPID_STR     ""
 #define ESCAPE_REVERSE_STR        ""
 #define ESCAPE_HIDDEN_STR         ""
 #define ESCAPE_CONCEAL_STR        ""

 #define ESCAPE_BRIGHT_OFF_STR     ""
 #define ESCAPE_BOLD_OFF_STR       ""
 #define ESCAPE_DIM_OFF_STR        ""
 #define ESCAPE_FAINT_OFF_STR      ""
 #define ESCAPE_ITALIC_OFF_STR     ""
 #define ESCAPE_UNDERLINE_OFF_STR  ""
 #define ESCAPE_BLINK_OFF_STR      ""
 #define ESCAPE_REVERSE_OFF_STR    ""
 #define ESCAPE_HIDDEN_OFF_STR     ""
 #define ESCAPE_CONCEAL_OFF_STR    ""

 #define ESCAPE_FGBLACK_STR        ""
 #define ESCAPE_FGRED_STR          ""
 #define ESCAPE_FGGREEN_STR        ""
 #define ESCAPE_FGYELLOW_STR       ""
 #define ESCAPE_FGBLUE_STR         ""
 #define ESCAPE_FGMAGENTA_STR      ""
 #define ESCAPE_FGCYAN_STR         ""
 #define ESCAPE_FGWHITE_STR        ""
 #define ESCAPE_FGGREY_STR        ""
 #define ESCAPE_FGLRED_STR       ""
 #define ESCAPE_FGLGREEN_STR       ""
 #define ESCAPE_FGLYELLOW_STR       ""
 #define ESCAPE_FGLBLUE_STR       ""
 #define ESCAPE_FGLMAGENTA_STR       ""
 #define ESCAPE_FGLCYAN_STR       ""
 #define ESCAPE_FGLWHITE_STR       ""

 #define ESCAPE_BGBLACK_STR        ""
 #define ESCAPE_BGRED_STR          ""
 #define ESCAPE_BGGREEN_STR        ""
 #define ESCAPE_BGYELLOW_STR       ""
 #define ESCAPE_BGBLUE_STR         ""
 #define ESCAPE_BGMAGENTA_STR      ""
 #define ESCAPE_BGCYAN_STR         ""
 #define ESCAPE_BGWHITE_STR        ""
 #define ESCAPE_BGGREY_STR       ""
 #define ESCAPE_BGLRED_STR       ""
 #define ESCAPE_BGLGREEN_STR       ""
 #define ESCAPE_BGLYELLOW_STR       ""
 #define ESCAPE_BGLBLUE_STR       ""
 #define ESCAPE_BGLMAGENTA_STR       ""
 #define ESCAPE_BGLCYAN_STR       ""
 #define ESCAPE_BGLWHITE_STR       ""

 #define FPRINTFTEXTALL(file,attr,fg,bg)  EMPTY
 #define FPRINTFTEXTATTR(file,attr)       EMPTY
 #define FPRINTFTEXTFG(file,fg)           EMPTY
 #define FPRINTFTEXTBG(file,bg)           EMPTY
 #define FPRINTFTEXTRESET(file)           EMPTY

 #define FPRINTFTEXT_ICH(file,n)          EMPTY
 #define FPRINTFTEXT_CUU(file,n)       EMPTY
 #define FPRINTFTEXT_CUD(file,n)          EMPTY
 #define FPRINTFTEXT_CUF(file,n)          EMPTY
 #define FPRINTFTEXT_CUB(file,n)          EMPTY
 #define FPRINTFTEXT_CNL(file,n)          EMPTY
 #define FPRINTFTEXT_CPL(file,n)          EMPTY
 #define FPRINTFTEXT_CHA(file,n)          EMPTY
 #define FPRINTFTEXT_CUP(file,n,m)        EMPTY
 #define FPRINTFTEXT_CHT(file,n)          EMPTY
 #define FPRINTFTEXT_ED(file,n)           EMPTY
 #define FPRINTFTEXT_EL(file,n)           EMPTY
 #define FPRINTFTEXT_IL(file,n)           EMPTY
 #define FPRINTFTEXT_DL(file,n)           EMPTY
 #define FPRINTFTEXT_DCH(file,n)          EMPTY
 #define FPRINTFTEXT_SU(file,n)           EMPTY
 #define FPRINTFTEXT_SD(file,n)           EMPTY
 #define FPRINTFTEXT_HVP(file,n,m)        EMPTY
 #define FPRINTFTEXT_SGR(file,n)          EMPTY
 #define FPRINTFTEXT_DSR(file)            EMPTY
 #define FPRINTFTEXT_SCP(file)            EMPTY
 #define FPRINTFTEXT_RCP(file)            EMPTY
 #define FPRINTFTEXT_SCA(file)            EMPTY
 #define FPRINTFTEXT_RCA(file)            EMPTY

 #define PRINTFTEXT_ICH(n)                EMPTY
 #define PRINTFTEXT_CUU(n)                EMPTY
 #define PRINTFTEXT_CUD(n)                EMPTY
 #define PRINTFTEXT_CUF(n)                EMPTY
 #define PRINTFTEXT_CUB(n)                EMPTY
 #define PRINTFTEXT_CNL(n)                EMPTY
 #define PRINTFTEXT_CPL(n)                EMPTY
 #define PRINTFTEXT_CHA(n)                EMPTY
 #define PRINTFTEXT_CUP(n,m)              EMPTY
 #define PRINTFTEXT_CHT(n)                EMPTY
 #define PRINTFTEXT_ED(n)                 EMPTY
 #define PRINTFTEXT_EL(n)                 EMPTY
 #define PRINTFTEXT_IL(n)                 EMPTY
 #define PRINTFTEXT_DL(n)                 EMPTY
 #define PRINTFTEXT_DCH(n)                EMPTY
 #define PRINTFTEXT_SU(n)                 EMPTY
 #define PRINTFTEXT_SD(n)                 EMPTY
 #define PRINTFTEXT_HVP(n,m)              EMPTY
 #define PRINTFTEXT_SGR(n)                EMPTY
 #define PRINTFTEXT_DSR()                 EMPTY
 #define PRINTFTEXT_SCP()                 EMPTY
 #define PRINTFTEXT_RCP()                 EMPTY
 #define PRINTFTEXT_SCA()                 EMPTY
 #define PRINTFTEXT_RCA()                 EMPTY

#endif /* _WITH_ANSI_ESCAPE_CODES */

#endif /* ESCAPE_CODES_H */
