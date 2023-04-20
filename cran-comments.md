R CMD gives no errors or warnings, and one note stating that this is a
new submission, which I cannot eliminate:

* checking CRAN incoming feasibility ... [8s] NOTE
Maintainer: 'William Lees <william@lees.org.uk>'

New submission

-----------------

I confirm that I have read and agree to the CRAN Repository Policy.

R CMD has beeen run on Windows, macOS and Linux using GitHub actions, and
also using win_devel.

-----------------

Thank you for earlier feedback. Since then I have:
- Proof-read the description text and made some changes
- added a reference describing the methods to the description text
- Replaced the use of T and F with TRUE and FALSE throughout the package
- Added examples to all exported functions
- Replaced print()/cat() with use message()/warning()  
- Removed all code that changed the user's current directory
- reduced the execution time of example code
- removed words that offended the spell checker

Best wishes
