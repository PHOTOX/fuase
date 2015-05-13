=================
How to contribute
=================

Discussion of ASE development takes place on the :ref:`ase-developer
<mailing_lists>` mailing list, and also sometimes on the #gpaw IRC
channel on freenode.

We welcome new developers who would like to help work on improving
ASE.  If you would like to contribute, your should first tell us what
you want to work on.  Use the mailing list for that.


SVN access
==========

We don't give new contributers write access to our SVN repository from
day one.  So, you will have to create a patch and send it to the
mailing list::

  $ svn checkout https://svn.fysik.dtu.dk/projects/ase/trunk myase
  $ cd myase
  $ # do your thing ...
  $ svn diff > patch.txt

Before you send the patch, *please* read our
:ref:`python_codingstandard` and learn how to use pep8.py,
pylint and epydoc:

* :ref:`pep8py`
* :ref:`pylint`
* :ref:`epydoc`

One of the current committers will look at the patch and give you some
feedback.  Maybe the patch is fine and the committer will commit it to
trunk.  There could also be some more work to do like:

* make it compatible with all supported pythons (see :ref:`download_and_install`).
* write more comments
* fix docstrings
* write a test
* add some documentation

Once everyone is happy, the patch can be applied.  This patch-feedback
loop is not something we have invented to prevent you from
contributing - it should be viewed as an opportunity for you to learn
how to write code that fits into the ASE codebase.  

After a couple of contributions, we will probably trust you enough to
add you as a committer.


Committers
==========

Here is the list of current committers:

==========  =========================  ===================================
user name   real name
==========  =========================  ===================================
anpet       Andrew Peterson            andrew_peterson:brown,edu
askhl       Ask Hjorth Larsen          asklarsen:gmail,com
bjork       Jonas Bjork                jonbj:ifm,liu,se
dlandis     David Landis               dlandis:fysik,dtu,dk
dulak       Marcin Dulak               dulak:fysik,dtu,dk
eojons      Elvar Örn Jónsson          elvar,jonsson:fysik,dtu,dk
ehermes     Eric Hermes                ehermes:chem,wisc,edu
gdonval     Gaël Donval                gael,donval:cnrs-imn,fr
getri       George Tritsaris           gtritsaris:seas,harvard,edu
grabow      Lars Grabow                grabow:uh,edu
hahansen    Heine Anton Hansen         hahansen:fysik,dtu,dk
ivca        Ivano Eligio Castelli      ivca:fysik,dtu,dk
jakobb      Jakob Blomquist            jakobb:fysik,dtu,dk
jber        Jon Bergmann Maronsson     jber:fysik,dtu,dk
jblomqvist  Janne Blomqvist            Janne,Blomqvist:tkk,fi
jensj       Jens Jørgen Mortensen      jensj:fysik,dtu,dk
jesperf     Jesper Friis               jesper,friis:sintef,no
jesswe      Jess Wellendorff Pedersen  jesswe:fysik,dtu,dk
jingzhe     Jingzhe Chen               jingzhe:fysik,dtu,dk
jkitchin    John Kitchin               jkitchin:andrew,cmu,edu
jussie      Jussi Enkovaara            jussi,enkovaara:csc,fi
kkaa        Kristen Kaasbjerg          kkaa:fysik,dtu,dk
kleis       Jesper Kleis               kleis:fysik,dtu,dk
kwj         Karsten Wedel Jacobsen     kwj:fysik,dtu,dk
lassebv     Lasse Vilhelmsen           lassebv:phys,au,dk
markus      Markus Kaukonen            markus,kaukonen:iki,fi
miwalter    Michael Walter             Michael,Walter:fmf,uni-freiburg,de
moses       Poul Georg Moses           poulgeorgmoses:gmail,com
mvanin      Marco Vanin                mvanin:fysik,dtu,dk
pastewka    Lars Pastewka              lars,pastewka:iwm,fraunhofer,de
s032082     Christian Glinsvad         s032082:fysik,dtu,dk
schenkst    Stephan Schenk             stephan,schenk:basf,com
schiotz     Jakob Schiotz              schiotz:fysik,dtu,dk
slabanja    Mattias Slabanja           slabanja:chalmers,se
strange     Mikkel Strange             mikkel,strange:gmail,com
tjiang      Tao Jiang                  tjiang:fysik,dtu,dk
tolsen      Thomas Olsen               tolsen:fysik,dtu,dk
==========  =========================  ===================================


Former committers:

==========  ======================  ===================================
anro        Anthony Goodrow         anro:fysik,dtu,dk 
carstenr    Carsten Rostgaard       carstenr:fysik,dtu,dk
hanke       Felix Hanke             F,Hanke:liverpool,ac,uk
s042606     Janosch Michael Rauba   s042606:fysik,dtu,dk
s052580     Troels Kofoed Jacobsen  s052580:fysik,dtu,dk
==========  ======================  ===================================
