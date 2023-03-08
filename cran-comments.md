R CMD gives no errors or warnings, and two notes:

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: 'William Lees <william@lees.org.uk>'
  
  New submission

❯ checking examples ... [36s] NOTE
  Examples with CPU (user + system) or elapsed time > 5s
                           user system elapsed
  generate_ogrdb_report   11.40   1.03   14.61
  genotype_statistics_cmd  4.68   0.69    6.29
  write_plot_file          4.29   0.56    6.16
  
I have made the example as minimal as I can.

Thank you for earlier feedback. I have:
- Proof-read the description text and made some changes
- added a reference describing the methods to the description text
- Replaced the use of T and F with TRUE and FALSE throughout the package
- Added examples to all exported functions
- Replaced print()/cat() with use message()/warning()  
- Removed all code that changed the user's current directory

Best wishes
