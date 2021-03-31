c          PostScript plot interface
      common/pltopn_info/ Munits,Lpopen,Ltrace,Margin
      integer Munits,Margin
      logical Lpopen,Ltrace

c * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - * - *
c Munits = currently selected ploter devices (not used)
c Lpopen = an open has been done (and data written to files, maybe)
c Ltrace = type subr. name and args at each entry point
c Margin = 0 plot full page
c          1 shrink plot to fit 1-inch margins on all four sides
c          >1 plotting routine should allow Margin inches top and right
