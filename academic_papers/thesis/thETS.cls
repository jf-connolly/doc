%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Modèle LaTeX pour le mémoire de thèse de l'ETS
%% Basé sur la classe de base book.cls, généralement positionné 
%% dans un système linux à
%% /usr/share/texmf/tex/latex/base/book.cls
%%
%% Modifications par Nicolas Gasset (nicolas.gasset.1@ens.etsmtl.ca)
%% [08 fev 2012]
%% - Tentative d amelioration des espaces dessus/dessous figure/tableau
%% - Commente \RequirePackage{mathptmx} (car incompatible avec formules mathematique en gras), et remplace par \RequirePackage{times}
%% - Supprime option "body={6.25in,8.75in}" lors de l appel du paquet geometry (c etait inutile/surimposé, cf. warning avec l option)
%%
%% Modifications par Nicolas Gasset (nicolas.gasset.1@ens.etsmtl.ca)
%% [03 fev 2012]
%% - Ajout du paquet hyperref
%% - Integration ecriture chapitre/annexe avec hyperref (cf modif Remi Menegon, Tennessee Carmel-Veilleux)
%% - Retablissement de la clef "hypertext" renommé a "hyperref"
%%
%% Modifications par Nicolas Gasset (nicolas.gasset.1@ens.etsmtl.ca)
%% [03 fev 2012]
%% -Ajout du compteur annexnum pour numerote les annexe en evitant d utiliser le compteur des chapitres
%% - Retablissemnt de la section "Numérotation des annexes en chiffres arabes" comme precedement
%%
%% Modifications par Nicolas Gasset (nicolas.gasset.1@ens.etsmtl.ca)
%% [03 fev 2012]
%% -Corrige bugs pour que le document soit en interline 1,5
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [25 Jan 2012]
%% Remplacement du logo Creative Commons en vectoriel, et déplacement des images
%% dans le dossier Figures.
%% (Merci à David Brassard pour les modifications)
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [23 Jan 2012]
%% Correction de "Directeur de thèse" pour les thèses par articles.
%% Ajout de la section "Revue de littérature pour les thèses par articles.
%% (Merci à Sébastien Rufiange pour les modifications)
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [30 Nov 2011]
%% Correction des marges pour correspondre aux normes de rédaction de l'ÉTS, en traffiquant avec fancyHDR
%% Correction de l'espacement float/texte en modifiant la variable \intextsep
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [11 Oct 2011]
%% Il n'est plus nécessaire d'utiliser les fonctions figureap et tableap en appendix.
%% (Merci à Felix Chenier pour la suggestion).
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [11 Oct 2011]
%% Simplification de la redéfinition des sous-sections des annexes:
%% Uniquement niveau 0 sur TOC, et pas d'affichage du N d'annexe dans les sous parties.
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [2 Sep 2011]
%% Nettoyage des différentes parties commentées pour plus de clareté.
%% Correction de la position du numéro de page selon les normes du gabarit. utilisation du package fancyhdr
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [31 Aout 2011]
%% Numérotation des figures: section.fig.
%% Suppression du ':' dans la caption des figures
%% Changement de la marge après retour à la ligne, pour les longues captions dans la liste de tables
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [13 Mars 2011]
%% Ligne 564: modification de l'affichage des titres des chapitres dans le TOC. Remplacement de la box par \numberline
%% Aussi, changement de la valeur de @tempdima pour que ça rende mieux.
%% suppression de la box avant le titre des tabs et figures dans les listes
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [12 Mai 2011]
%% \setcounter{page}{1} dans l'intro pour qu'elle soit en page 1
%% Suppression du reset de page dans title et jury pour que le titre soit la page I
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [11 Mai 2011]
%% Utilisation d'une newcommand plutôt qu'une savebox pour stocker les mots clé des abstract.
%% Corrige le problème du non retour à la ligne.
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [30 Mars 2011]
%% modification de l'affichage des titres des chapitres dans le TOC. Le nom et numéro sont maintenant dans 
%% une box de 7em, ce qui aligne les titres. 
%% changement de la notation des figures: Exemple Figure-A V-1. Aussi, pas en gras.
%% Changement de l'ordre d'apparition du jury: examinateur externe avant membre du jury.
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [24 Mars 2011]
%% Ajout de ':' après 'jury composé de'
%% Ajout de 3 retours à la ligne avant l'affichage des mots clés pour les abstract.
%% Commentaire de la section pour restaurer l'indentation des puces
%% ``Figure'' plus en gras.
%%  Ajout de saut de ligne avant figure.
%%
%% Modifications par Christophe Pagano (cpagano@livia.etsmtl.ca)
%% [9 Mars 2011]
%% Changement de la taille du texte (small) et de l'interligne (simple) pour la partie creative commons et copyright.
%%
%% Modifications par Paulo Cavalin (cavalin@livia.etsmtl.ca)
%% [17 dec 2010]
%% Modifications pour l'impression recto-verso
%%
%% [16 nov 2010]
%% Traduction en anglais
%%
%% [27 aout 2010]
%% Mises à jour du gabarit pour qu'il soit disponibilisé officielement par le décanat
%% Utilisation du package geometry pour définir les marges
%% corrections relatives au guide daté du 27-08-2010
%%
%% Corrections apportées par Cécile Le Cocq (cecile.lecocq@etsmtl.ca)
%%
%% [17 mars 2008]
%% corrections relatives au nouveau guide de rédaction daté du 12 novembre 2007
%%
%% Commentaires (et modifications?) par Luis Da Costa (costa@livia.etsmtl.ca)
%%
%% [14 Juillet 2006]  
%% Il faut laisser une ligne en blanc AVANT et APRES une citation de plus de 3 lignes: 
%% 
%% Zinedine Zidane a remarque 
%% [LIGNE EN BLANC]
%% \begin{quote} 
%%      parler de ma mere et de ma soeur merite bien un coup de boule. 
%% \end{quote} 
%% [LIGNE EN BLANC]
%%  
%% [14 Juillet 2006] 
%% Dans la table des matieres il doit y avoir une indentation. 
%% J'ai donc commenté ce que Seb avait implementé (chercher LUIS dans ce document
%% pour voir).  
%% 
%% Version 0.4 
%% Date : 2005/08/04
%%
%% Modifications par S. Piccand (sebastien.piccand@livia.etsmtl.ca)
%%
%% thETS4 est basé sur thETS3
%% Quelques corrections au niveau de la table des matières
%%  - la division en section d'un chapitre est séparée par une ligne de l'intitulé du chapitre
%%  - les annexes sont correctement gérées
%%  - la bibliographie est corrigée en utilisant le package natbib.sty et le style ETS.bst
%%  - ajout de l'option memoire. Par defaut il s'agit d'une thèse
%%  - correction des figures, sous figures, tables
%%  - correction des footnotes, utiliser \footnote{} comme d'habitude
%%  - correction de l'environnement "enumerate"
%%  - correction des citations utiliser \quote{} comme d'habitude
%%
%% Macros supplémentaires pour ajouter plus facilement des figures/tableaux et y faire référence
%% 
%% Non corrigé de la version 0.3 :
%%  * Dans les différentes listes (figures, tables, etc.) il n'y a pas assez d'espace entre la fin
%%     du titre et les numéros. Une solution consiste à utiliser des noms de figure courts.
%%  * Les titres des sections/chapitre/soussections/etc doivent être courts pour ne pas passer à
%%     la ligne dans la table des matières.
%%  * Les citations ne sont pas encore bonnes.
%%
%%
%% Version 0.3
%% Date : 2004/12/17
%% 
%% Modifications par G. Tremblay (guillaume.tremblay@livia.etsmtl.ca)
%% 
%% thETS3 est basé sur thETS2 
%% Quelques corrections apportées par rapport à la version 0.2: 
%% 
%% Utilisation de algorithme2e.sty sans problèmes de compatibilité. Correction faite par F. Grandidier
%% et appliquée par G. Tremblay
%% 
%% Les citations de plus de 3 lignes sont à interligne simple tel qu'exigé par le Décanat (utilisez begin{quoteETS})
%%
%% Les footnotes sont remplacées par de footnoteETS pour ne pas avoir d'indentation de paragraphe
%%
%%
%% Version 0.2 
%% Date: 2003/05/06
%% 
%% Par F. Grandidier (gfrede@livia.etsmtl.ca)
%%
%% Ajout ver 0.2: alignement du titre du chapitre dans la table
%%                des matières lorsqu'il tient sur plusieurs
%%                lignes
%% 
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\filedate{2003/01/24}

\NeedsTeXFormat{LaTeX2e}[1994/06/01] %% Mis par défaut
\ProvidesClass{thETS}[\filedate\space Standard ETS class]

%% Mêmes options possibles que celle de la classe book.cls
%% Cependant tous ne va peut être pas marcher comme dans book...
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}} 

\newcommand{\TypeDoc}{empty}

\DeclareOption{these}{%
	\renewcommand{\TypeDoc}{these}
}
	
\DeclareOption{memoire}{%
  \renewcommand{\TypeDoc}{memoire}
}
\DeclareOption{thesepararticles}{%
  \renewcommand{\TypeDoc}{thesepararticles}
}

\newcommand{\TypeLic}{reserved}
\DeclareOption{creativecommons}{
   \renewcommand{\TypeLic}{creativecommons}
}

%% Ajout de l'option web pour activer hyperref (Menegon 26-07-2010)
\newcommand{\WithHyperref}{false}
\DeclareOption{hyperref}{%
  \renewcommand{\WithHyperref}{true}
}

%% Margin variables for fancy HDR

\newcommand{\fanctmargin}{0.39in}
\newcommand{\fanchsep}{0.56in}
\newcommand{\newbmargin}{0.58in}  %(0.75 - \fanchsep - \fanctmargin)

\ProcessOptions

\LoadClass{book} %% Chargement de la classe book
%\RequirePackage{mathptmx} %% Police Times .... mais incompatible avec lettre en gras dans equations! utilisation de times juste?!
\RequirePackage{times}

%\RequirePackage[frenchb]{babel} %% Francisation du document
\RequirePackage[english]{babel} %% Francisation du document
\RequirePackage[utf8]{inputenc} %% Prise en compte de la frappe de caractères accentués
\RequirePackage[T1]{fontenc} %% Codage standard des fontes LaTeX pour le français
\RequirePackage{calc} %% permet de réaliser des additions et soustractions de valeurs lors d'affectation de variables
\RequirePackage{setspace}
\RequirePackage{ifthen} %% Permet de réaliser des tests sur des variables
\RequirePackage{amsthm} % clc
\RequirePackage{longtable}
\RequirePackage{float}
\RequirePackage{listings} %% Code listings (extraits de code)
\RequirePackage[labelfont={bf},margin=10pt,justification=centering]{caption}
\RequirePackage[paper=letterpaper,tmargin=1.5in,lmargin=1.5in,rmargin=0.75in,bmargin=\newbmargin]{geometry}
\RequirePackage{fancyhdr}
\RequirePackage[caption=false]{subfig} %% Utilisation des subfigures
\RequirePackage{tabto}
\setlength{\headheight}{15pt} % pour supprimer un warning avec fancyhdr

%% Exec de l'option web pour activer hyperref (Menegon 26-07-2010)
\ifthenelse{\equal{\WithHyperref}{true}}{
  \RequirePackage[%
pdfborder={0 0 0 0},%
%%%%               Bookmarck dans le pdf
bookmarks=true,%                  %% Cree les bookmarks ...
bookmarksnumbered=false,%          %% ... avec des numro ...
bookmarksopen=true,%              %% ... et les affiche a l'ouverture du pdf
bookmarksopenlevel=2,%            %%     jusqu au niveau
%%%%               Options pour le rendu
pdfhighlight=/O%                 %% souligne les liens lors d'un clic (defaut /I : inverse)
]{hyperref}
}{
\RequirePackage[hyphens]{url}
\newcommand{\href}[2]{\url{#1}}
\newcommand{\texorpdfstring}[2]{#1}
}
\urlstyle{rm} % Hyphenation des références

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Définition de commandes pour la présentation du jury
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\PrefDirect}{Mr.}
\newcommand{\NomDirect}{Nom du directeur de recherche}
\ifthenelse{\equal{\TypeDoc}{these}}
{\newcommand{\TxtDirect}{Thesis Director}}
{\ifthenelse{\equal{\TypeDoc}{thesepararticles}}
	{\newcommand{\TxtDirect}{Thesis Director}}
	{\newcommand{\TxtDirect}{memorandum director}}
}
\newcommand{\AffilDirect}{Affiliation du directeur}
\newcommand{\PrefCoDirect}{Mr.}
\newcommand{\NomCoDirect}{Nom du co-directeur}
\newcommand{\TxtCoDirect}{Thesis co-Director}
\newcommand{\AffilCoDirect}{Affiliation du codirecteur}
\newcommand{\PrefPresid}{Mr.}
\newcommand{\NomPresid}{Nom du pr\'esident du jury}
\newcommand{\TxtPresid}{Committee President}
\newcommand{\AffilPresid}{Affiliation du pr\'esident du jury}
\newcommand{\PrefExamExt}{Mr.}
\newcommand{\NomExamExt}{Nom de l'examinateur externe}
\newcommand{\TxtExamExt}{External Examiner}
\newcommand{\PlusTxtExamExt}{Blabla}
\newcommand{\AffilExamExt}{Affiliation de l'examinateur externe}
\newcommand{\PrefJury}{Mr.}
\newcommand{\NomJury}{Nom d'un membre du jury}
\newcommand{\TxtJury}{Examiner}
\newcommand{\PlusTxtJury}{Blabla}
\newcommand{\AffilJury}{Affiliation du membre du jury}
\newcommand{\PrefJuryB}{Mr.}
\newcommand{\NomJuryB}{Nom d'un membre du jury}
\newcommand{\TxtJuryB}{invited examiner}
\newcommand{\PlusTxtJuryB}{Blabla}
\newcommand{\AffilJuryB}{Affiliation du membre du jury}
  
%% Directeur de thèse
\newcommand{\directeur}[3]{%
  \renewcommand{\PrefDirect}{#1}%
  \ifthenelse{\equal{#1}{Mrs.}}{% Permet de modifier pour le féminin
    \ifthenelse{\equal{\TypeDoc}{these}}
    {\renewcommand{\TxtDirect}{Thesis Director}}
		{\ifthenelse{\equal{\TypeDoc}{thesepararticles}}
			{\renewcommand{\TxtDirect}{Thesis Director}}
			{\renewcommand{\TxtDirect}{memorandum director}}}}
	{}%
  \renewcommand{\NomDirect}{#2}%
  \renewcommand{\AffilDirect}{#3}%
  }%  
%%% Co-directeur de thèse
\newcommand{\codirecteur}[3]{%
  \renewcommand{\PrefCoDirect}{#1}%
  \ifthenelse{\equal{#1}{Mme}}{%
    \renewcommand{\TxtCoDirect}{codirector}}{}%
  \renewcommand{\NomCoDirect}{#2}%
  \renewcommand{\AffilCoDirect}{#3}%
  }%
%% Président du jury
\newcommand{\president}[3]{%
  \renewcommand{\PrefPresid}{#1}%
  \ifthenelse{\equal{#1}{Mme}}{%
    \renewcommand{\TxtPresid}{committee president}}{}%
  \renewcommand{\NomPresid}{#2}%
  \renewcommand{\AffilPresid}{#3}%
  }%
%% Examinateur externe
\newcommand{\examinexterne}[4]{%
  \renewcommand{\PrefExamExt}{#1}%
  \ifthenelse{\equal{#1}{Mme}}{%
    \renewcommand{\TxtExamExt}{extern examiner}}{}%
  \renewcommand{\NomExamExt}{#2}%
  \renewcommand{\AffilExamExt}{#3}%
  \renewcommand{\PlusTxtExamExt}{#4}%
  }%
%% Membre du jury
\newcommand{\jury}[4]{%
  \renewcommand{\PrefJury}{#1}%
  \ifthenelse{\equal{#1}{Mme}}{%
    \renewcommand{\TxtJury}{examiner}}{}%
  \renewcommand{\NomJury}{#2}%
  \renewcommand{\AffilJury}{#3}%
  \renewcommand{\PlusTxtJury}{#4}%
  }%
%% Membre du jury
\newcommand{\juryB}[4]{%
  \renewcommand{\PrefJuryB}{#1}%
  \ifthenelse{\equal{#1}{Mme}}{%
    \renewcommand{\TxtJuryB}{invited examiner}}{}%
  \renewcommand{\NomJuryB}{#2}%
  \renewcommand{\AffilJuryB}{#3}%
  \renewcommand{\PlusTxtJuryB}{#4}%
  }%
%% Bon là, si il y a d'autres membres ça devient le bordel. 2 possibilité : 
%% - rendre le code dynamique, mais ça à l'air vachement compliqué
%% - créer des variables examext2, juryB, jury3,... (plus facile)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pour le copyright, c'est Prénom Nom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\authorcopyright#1{\gdef\@authorcopyright{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable associée aux dates de
%% soutenance et de dépôt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\datesoutenance#1{\gdef\@datesoutenance{#1}}
\def\datedepot#1{\gdef\@datedepot{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable pour la mise en page des chapitres dans
%% la table des matières
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newlength{\largChap}
\newlength{\largNum}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Marges autour des objets flottants : figures et tableaux
%% Ce n'est pas au point, il y aurait des ajustements à faire.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setlength\abovecaptionskip{\baselineskip}
\setlength\belowcaptionskip{\baselineskip}
\setlength\textfloatsep{\baselineskip}
\setlength\intextsep{\baselineskip}

%%%%%%%%%%%%%%%%%%%
%% Marge la table des matières et les différentes listes
%%%%%%%%%%%%%%%%%%%
\newlength{\tocRightMargin} % Marge droite
\setlength{\tocRightMargin}{2cm}
\newlength{\tocLeftMarginSecondLine}
\setlength{\tocLeftMarginSecondLine}{6.1em} % Marge gauche quand on va à la ligne pour les sections et sous sections
\newlength{\tocLeftMarginSecondLineChapter}
\setlength{\tocLeftMarginSecondLineChapter}{6.1em} % Marge gauche quand on va à la ligne pour les chapitres et annexes
\newlength{\tocLeftMarginSecondLineFigure}
\setlength{\tocLeftMarginSecondLineFigure}{6.5em} % Marge gauche quand on va à la ligne pour les figures
\newlength{\tocLeftMarginSecondLineTable}
\setlength{\tocLeftMarginSecondLineTable}{6em} % Marge gauche quand on va à la ligne pour les tableaux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fraction maxi des objets flottants 
%% en haut et bas de page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\floatpagefraction}{.99}
\renewcommand{\textfraction}{.01}
\renewcommand{\topfraction}{.99}
\renewcommand{\bottomfraction}{.99}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interligne du document
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\mystretch{1.5}
\def\baselinestretch{\mystretch}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Espace entre 2 paragraphes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setlength{\parskip}{12pt} % 12pt 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pas d'indentation des paragraphes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setlength{\parindent}{0pt}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Pour éviter les veuves et orphelins
%% (lignes seules en début ou fin de pages)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\widowpenalty=10000
\clubpenalty=10000
\raggedbottom

%%% footnote 
%
% Pour enlever l'indentation des footnotes (LUIS)
%
\newlength{\NameOfNewDefinedLength}

\newcommand{\footnoteETS}[1]%
{\footnote{\setlength{\NameOfNewDefinedLength}{\textwidth}
                 \advance\leftskip -2.7mm  %%
                 \addtolength{\NameOfNewDefinedLength}{-2em+2.mm}
                 \hskip -1.mm
                 \parbox[t]{\NameOfNewDefinedLength}{\raggedright #1}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ajout d'un compteur pour les annexes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcounter{annexnum}           
\setcounter{annexnum}{-1}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Page titre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\lediplome}{
  THE DEGREE OF DOCTOR OF PHILOSOPHY\\
  Ph.D.
}
\newcommand{\diplome}[1]{
  \renewcommand{\lediplome}{#1}
}  

\renewcommand{\maketitle}{\begin{titlepage}
 \renewcommand{\baselinestretch}{1}\small\normalsize % Interligne simple
    \let\footnotesize\small%
    \let\footnoterule\relax%
    \let \footnote \thanks%
   \begin{center}
      ÉCOLE DE TECHNOLOGIE SUPÉRIEURE\\
      UNIVERSITÉ DU QUÉBEC
      \vfill
      \ifthenelse{\equal{\TypeDoc}{thesepararticles}}{
        MANUSCRIPT-BASED THESIS PRESENTED
      }{
        THESIS PRESENTED
      } TO\\
      ÉCOLE DE TECHNOLOGIE SUPÉRIEURE
      \vfill
      IN PARTIAL FULFILLMENT OF THE REQUIREMENTS FOR\\
      \lediplome 
      \vfill
      BY\\
	\@author
      \vfill
      \MakeUppercase{\@title}%
      \vfill
      MONTR\'{E}AL, \MakeUppercase{\@datedepot}%
      
      \vspace{1em}
      {\scriptsize 
      \ifthenelse{\equal{\TypeLic}{reserved}}{ 
          \copyright~Copyright \number\year \ reserved by \@authorcopyright
      }{
         \includegraphics[width=50bp]{Figures/logocc} \@authorcopyright \ \number\year
      }
      }%
   \end{center}
   \end{titlepage}
   \clearpage

   %Verso 
  \newpage     
  ~\\
  \vspace{15cm}
   
   \ifthenelse{\equal{\TypeLic}{reserved}}{%
   \fbox{\parbox{6.125in}{ %% Anciennement 6.4
    %% Interligne de 1 et taille de 10pt pour le texte
   \renewcommand{\baselinestretch}{1}\footnotesize\footnotesize

   \copyright~Copyright reserved
    \newline
   It is forbidden to reproduce, save or share the content of this document either in whole or in parts. The reader who 
   wishes to print or save this document on any media must first get the permission of the author.
   }}
   }{%
   \fbox{\parbox{6.125in}{
    %% Interligne de 1 et taille de 10pt pour le texte
   \renewcommand{\baselinestretch}{1}\footnotesize\footnotesize

   \includegraphics[width=50bp]{Figures/logocc}
   \newline
   This Creative Commons license allows readers to download this work and share it with others as long as the author is 
   credited. The content of this work cannot be modified in any way or used commercially.
   }}
   }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Page de présentation du jury
% LUIS A ELIMINE L'ESPACE ENTRE LES NOMS ET LES AFFILIATIONS (Décanat...) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand\presentjury{\begin{jurypage}%
  \let\footnotesize\small%
  \let\footnoterule\relax%
  \begin{center}%
    \renewcommand{\baselinestretch}{2}\small\normalsize
    \textbf{BOARD OF EXAMINERS}\\
    \ifthenelse{\equal{\TypeDoc}{these}}{
      THIS THESIS HAS BEEN EVALUATED
    }{THIS THESIS HAS BEEN EVALUATED
    }\\
    BY THE FOLLOWING BOARD OF EXAMINERS:\\
    \vspace{6em}%
  \end{center}
  \begin{flushleft}
    \renewcommand{\baselinestretch}{1}\small\normalsize
    % Directeur de mémoire
    \PrefDirect~\NomDirect,~\TxtDirect\\
    \AffilDirect\par
    \vspace{1em}
    % Codirecteur optionnel
    \ifthenelse{\equal{\NomCoDirect}{Nom du co-directeur}}{}{%
      \PrefCoDirect~\NomCoDirect,~\TxtCoDirect\\%
      \vspace{3pt}
      \AffilCoDirect\par%
      \vspace{1em}
    }%
    % Président du jury
    \PrefPresid~\NomPresid,~\TxtPresid\\%
    \AffilPresid\par%
    \vspace{1em}%

    % Examinateur externe
    \ifthenelse{\equal{\NomExamExt}{Nom de l'examinateur externe}}{}{%
      \PrefExamExt~\NomExamExt,~\TxtExamExt\\
      \AffilExamExt\par%
      \vspace{1em}
    }%
    
    
 % Jury
    \ifthenelse{\equal{\NomExamExt}{Nom d'un membre du jury}}{}{%
	\PrefJury~\NomJury,~\TxtJury\\
      \AffilJury\par%
      \vspace{1em}
    }%


   
    
    % JuryB optionnel
    \ifthenelse{\equal{\NomJuryB}{Nom d'un membre du jury}}{}{%
      \PrefJuryB~\NomJuryB\ifthenelse{\equal{\PlusTxtJuryB}{}}{,~}{,~\PlusTxtJuryB}%
      \TxtJuryB\\%
      \AffilJuryB\par%
      \vspace{1em}
    }%
      
  \end{flushleft}
  \begin{center}%
    \vspace{3em}%
    \renewcommand{\baselinestretch}{2}\small\normalsize%
    \ifthenelse{\equal{\TypeDoc}{these}}{
      THIS THESIS
    }{
      THIS THESIS
    } WAS PRESENTED AND DEFENDED\\
      IN THE PRESENCE OF A BOARD OF EXAMINERS AND PUBLIC\\
     ON \MakeUppercase{\@datesoutenance}\\
    AT ÉCOLE DE TECHNOLOGIE SUPÉRIEURE
  \end{center}\par%
  \vfil\null%
  \end{jurypage}%
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Profondeur de la numérotation des divisions des chapitres
%% Valeur de 1(\section) à 5(\subparagraph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setcounter{secnumdepth}{5}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Changement dans la table des matières 
%% des entrées associées aux chapitres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\def\@chapter[#1]#2{
  \ifnum \value{annexnum}=-1%
     \ifnum \c@secnumdepth >\m@ne%
       \if@mainmatter%
         \refstepcounter{chapter}%
         \typeout{\@chapapp\space\thechapter}%
         %\addcontentsline{toc}{chapter}{\protect\makebox[7em][l]{\protect\MakeUppercase{\@chapapp\ \thechapter}} \protect\MakeUppercase{#1} }%
         %      \addcontentsline{toc}{chapter}{\protect\numberline{\protect\MakeUppercase{\@chapapp\ \thechapter}}\protect\MakeUppercase{#1}}%
         
         % Tennessee Carmel-Veilleux (25 juillet 2010): Modif pour viter les erreurs de \uppercase avec hyperref
%         \@ifundefined{texorpdfstring}{
         \ifthenelse{\equal{\WithHyperref}{false}}{
            % Hyperref non dtect: on fait comme si de rien tait
   %         \addcontentsline{toc}{chapter}{\protect{\uppercase{\@chapapp\ \thechapter\hskip 1em #1}}}%
            \addcontentsline{toc}{chapter}{\protect\numberline{\protect\MakeUppercase{\@chapapp\ \thechapter}}\protect\MakeUppercase{#1}}%
        
         }{
            % Hyperref dtect: deux versions du label: une premire formate pour TeX, une deuxime simplifie pour les bookmarks
   %         \addcontentsline{toc}{chapter}{\protect{\texorpdfstring{\uppercase{\@chapapp\ \thechapter\hskip 1em #1}}{\@chapapp\ \thechapter\ #1}}}%
            \addcontentsline{toc}{chapter}{\protect{\texorpdfstring{\numberline{\protect\MakeUppercase{\@chapapp\ \thechapter}}\protect\MakeUppercase{#1}}{\@chapapp\ \thechapter\ #1}}}%
         }
       \else%
         \addcontentsline{toc}{chapter}{#1}%
       \fi%
     \else%
       \addcontentsline{toc}{chapter}{#1}%
     \fi%
  \else%     
     \ifnum \c@secnumdepth >\m@ne%
       \if@mainmatter%
         \refstepcounter{annexnum}%
         \typeout{\@chapapp\space\thechapter}%
%         \@ifundefined{texorpdfstring}{
         \ifthenelse{\equal{\WithHyperref}{false}}{
          % Hyperref non détecté: on fait comme si de rien était
      %      \addcontentsline{toc}{chapter}{\protect{\uppercase{\@chapapp\ \thechapter\hskip 1em #1}}}%
            \addcontentsline{toc}{chapter}{\protect\numberline{\protect\MakeUppercase{\@chapapp\ \thechapter}}\protect\MakeUppercase{#1}}%
         }{
          % Hyperref détecté: deux versions du label: une première formatée pour TeX, une deuxième simplifiée pour les bookmarks
      %     \addcontentsline{toc}{chapter}{\protect{\texorpdfstring{\uppercase{\@chapapp\ \thechapter\hskip 1em #1}}{\@chapapp\ \thechapter\ #1}}}%
           \addcontentsline{toc}{chapter}{\protect{\texorpdfstring{\numberline{\protect\MakeUppercase{\@chapapp\ \thechapter}}\protect\MakeUppercase{#1}}{\@chapapp\ \thechapter\ #1}}}%
         }
     \else%
         \addcontentsline{toc}{chapter}{#1}%
       \fi%
     \else%
       \addcontentsline{toc}{chapter}{#1}%
     \fi%
  \fi%
  \chaptermark{#1}%
  \@makechapterhead{#2}%
  \@afterheading%
  }%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modification des chapitres numérotés
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\@makechapterhead}[1]{{
        \normalfont\bfseries%
        \begin{center}
        \etchapterheadstyle{\@chapapp{ \thechapter}} % CHAPTER 1
        \vskip 3\p@%
        \etchapterheadstyle{#1} % TITLE OF CHAPTER
        \end{center}
        \vskip 5\p@%
        \thispagestyle{empty}%
}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modification des chapitres sans numéro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\@makeschapterhead}[1]{{
    \thispagestyle{empty}
        \normalfont\bfseries%
        \begin{center}
        \etchapterheadstyle{#1} % TITLE OF CHAPTER
        \end{center}
}}

\newcommand{\etchapterheadstyle}[1]{\MakeUppercase{#1}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Définition des espacements avant et après les différentes
%% divisions de chapitre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\section{\@startsection{section}{1}{\z@}%
                                   {-0.25ex \@plus -1ex \@minus -.2ex}%
                                   {0.5ex \@plus.2ex}%
                                   {\normalfont\normalsize\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-0.25ex\@plus -1ex \@minus -.2ex}%
                                     {0.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-0.25ex\@plus -1ex \@minus -.2ex}%
                                     {0.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {-0.25ex \@plus -1ex \@minus -.2ex}%
                                    {0.5ex \@plus .2ex}%
                                    {\normalfont\normalsize\bfseries}}
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {0.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\normalfont\normalsize\bfseries}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rédéfinition de l'environnement de la page de titre 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewenvironment{titlepage}{%
  \if@twocolumn
     \@restonecoltrue\onecolumn%
  \else
     \@restonecolfalse\newpage
  \fi
  \thispagestyle{empty}%
  \setcounter{page}{1}
  }%
{\if@restonecol\twocolumn \else \newpage \fi%
  \if@twoside\else%
  \setcounter{page}\@ne%
  \fi
  \thispagestyle{empty}%
  }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Définition de l'environnement de la présentation du jury
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{jurypage}{%
  \if@twocolumn
     \@restonecoltrue\onecolumn%
  \else
     \@restonecolfalse\newpage
  \fi
  \thispagestyle{empty}%
  }%
{\if@restonecol\twocolumn \else \newpage \fi%
  \if@twoside\else%
  \setcounter{page}\@ne%
  \fi
  \thispagestyle{empty}%
  \cleardoublepage
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Abstract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{abstract}[1]{%
  \newcommand\MotCles{#1}
  \if@twocolumn
     \@restonecoltrue\onecolumn%
  \else
     \@restonecolfalse\newpage
  \fi
  \cleardoublepage
  \thispagestyle{empty}% Pas de numérotation sur cette page
  \renewcommand{\baselinestretch}{1}\small\normalsize % Interligne simple
  \begin{center}
    \normalfont%
    \MakeUppercase{\bfseries\@title}%
    \vskip 6\p@%
    \@author%
    \vskip 6\p@%
    \MakeUppercase{\textbf{ABSTRACT}}%
  \end{center}
  
}
{
  \vspace{2\baselineskip} %Ajout de 3 retours à la ligne avant les Mot clés

  {\bf Keywords:} \MotCles

  \if@restonecol\twocolumn \else \newpage \fi%
  \cleardoublepage
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sommaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{sommaire}[2]{%
  \newcommand\KeyWords{#2}
  \if@twocolumn
     \@restonecoltrue\onecolumn%
  \else
     \@restonecolfalse\newpage
  \fi
  \cleardoublepage
  \thispagestyle{empty}% Pas de numérotation sur cette page
  \renewcommand{\baselinestretch}{1}\small\normalsize % Interligne simple
  \begin{center}
    \normalfont%
    \MakeUppercase{\bfseries#1}%
    \vskip 6\p@%
    \@author%
    \vskip 6\p@%
    \MakeUppercase{\textbf{RÉSUMÉ}}%
  \end{center}
}
{
  \vspace{3\baselineskip} %Ajout de 3 retours à la ligne avant les Mot clés
  {\bf Mots-clés:} \KeyWords

  \if@restonecol\twocolumn \else \newpage \fi%
  \cleardoublepage
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Avant-propos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{avantpropos}{%
  \chapter*{Avant-propos}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Introduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{introduction}{%
  \cleardoublepage
  \chapter*{Introduction}
  \addcontentsline{toc}{chapter}{INTRODUCTION} % Ajout dans la table des matières
  \setcounter{page}{1}
 }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Revue de littérature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{revuedelitterature}{%
  \cleardoublepage
  \chapter*{Literature review}
  \addcontentsline{toc}{chapter}{LITERATURE REVIEW} % Ajout dans la table des matières
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conclusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{conclusion}{%
  \chapter*{Conclusion}
  \addcontentsline{toc}{chapter}{CONCLUSION} % Ajout dans la table des matières
 }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remerciement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{acknowledgements}{%
  \chapter*{Acknowledgements}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Commande pour écrire ``Page'' 
%% dans les différentes tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\tdmpage}[0]{\hspace{\fill}\mbox{Page}\vskip 0em\par}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table des matières
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\renewcommand\tableofcontents{%
%  \chapter*{\contentsname}
%  \tdmpage
%  {
%    % Espace avant les différentes lignes de la table des matières
%    \setlength{\parskip}{0pt}%
%    % Interligne simple pour la table des matières
%    \renewcommand{\baselinestretch}{1}\small\normalsize%
%    \@starttoc{toc}
%    }
%}

%%% Modification de l'espace entre les points de la table des matières
\renewcommand{\@dotsep}{1} % (\@dotsep est une valeur)


% Les parties n'étant pas utilisées on s'en sert pour pouvoir ajouter annexe correctement dans la table des matieres
\renewcommand*\l@part[2]{%
  \settowidth{\largChap}{\MakeUppercase{\@chapapp\ \thechapter\hskip 1em}}
  \ifnum \c@tocdepth >\m@ne%
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@%
    \setlength\@tempdima{\largChap}%
    \begingroup
      \parindent \z@ \rightskip \tocRightMargin%
      \parfillskip -\tocRightMargin%
      \leavevmode %
      \advance\leftskip \tocLeftMarginSecondLineChapter%
      \hskip -\leftskip%
      #1\nobreak\leaders\hbox{$\m@th\mkern \@dotsep mu\hbox{}\mkern \@dotsep mu$}\hfil \nobreak\hb@xt@\@pnumwidth{\hss }\par%
      \penalty\@highpenalty%
    \endgroup
  \fi
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modif pour ne pas avoir le chapitre en gras 
%% dans la table des matières
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand*\l@chapter[2]{%
  \settowidth{\largChap}{\MakeUppercase{\@chapapp\ \thechapter\hskip 1em}}
  \settowidth{\largNum}{\hss #2}
  \ifnum \c@tocdepth >\m@ne%
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@%
    \setlength\@tempdima{\largChap}%
    \begingroup
      \parindent \z@ \rightskip \tocRightMargin%
      \parfillskip -\tocRightMargin%
      \leavevmode %
      \advance\leftskip \tocLeftMarginSecondLine%
      \hskip -\leftskip%
      #1\nobreak\leaders\hbox{$\m@th \hbox{.}\mkern \@dotsep mu$}\hfil \nobreak\hb@xt@%
      \largNum{\hss #2}\par%
      \penalty\@highpenalty%
    \endgroup
  \fi}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table des matières
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\tableofcontents{%
%  \renewcommand{\baselinestretch}{1}\small\normalsize%
  \chapter*{\contentsname}
  %\vskip 1\baselineskip
  \tdmpage
    {
      % Espace avant les différentes lignes de la table des matières
      \setlength{\parskip}{0pt}%
      % Interligne simple pour la table des matières
      \renewcommand{\baselinestretch}{1}\small\normalsize%
      \@starttoc{toc}
    }
}
\renewcommand*\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne%
    \addpenalty{-\@highpenalty}%
    \vskip 1em \@plus\p@%
    \setlength\@tempdima{7.0em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth%
      \parfillskip -\@pnumwidth%
      \leavevmode %
      \advance\leftskip\@tempdima%
      \hskip -\leftskip%
      #1\nobreak
      \leaders\hbox{$\m@th\mkern \@dotsep mu\hbox{.}\mkern \@dotsep mu$}
      \hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par%
      \penalty\@highpenalty%
    \endgroup
  \fi
  }
\renewcommand*\l@section{\@dottedtocline{1}{0em}{3.5em}}
\renewcommand*\l@subsection{\@dottedtocline{2}{3.5em}{3.5em}}
\renewcommand*\l@subsubsection{\@dottedtocline{3}{7em}{3.5em}}
\renewcommand*\l@paragraph{\@dottedtocline{4}{2.5em}{3.5em}}
\renewcommand*\l@subparagraph{\@dottedtocline{5}{2.5em}{3.5em}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Profondeur des entrées dans la table 
%% des matières: 3=\subsubsection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\setcounter{tocdepth}{3}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liste des figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\listoffigures{%
%    \renewcommand{\baselinestretch}{1}\small\normalsize%
    \chapter*{\listfigurename}
    %\vskip 1\baselineskip
    \tdmpage
    {
      \setlength{\parskip}{0pt}%
      % Interligne simple pour la table des matières
      \renewcommand{\baselinestretch}{1}\small\normalsize%
      \@starttoc{lof}%
    }
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numérotation des figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand \thefigure
     {\arabic{chapter}.\@arabic\c@figure}
\renewcommand*\l@figure[2]{%
  \settowidth{\largNum}{\hss #2}
  \ifnum \c@tocdepth >\m@ne%
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@%
    \setlength\@tempdima{3.4em}%
        \noindent%
    \begingroup
      \pretolerance=10000
      \parindent \z@ \rightskip \tocRightMargin%
      \parfillskip -\tocRightMargin%
     \leavevmode \normalsize%
     \advance\leftskip \tocLeftMarginSecondLineFigure%
      \hskip -\leftskip%
      {\figurename\mbox{\hspace{4pt}}#1}\nobreak%
       \leaders\hbox{$\m@th%
        \mkern \@dotsep mu\hbox{.}\mkern \@dotsep%
        mu$} \hfil\nobreak\hb@xt@%
                    \largNum{\hss #2}\par%
      \penalty\@highpenalty%
    \endgroup
  \fi}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liste des tableaux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\listoftables{%
%  \renewcommand{\baselinestretch}{1}\small\normalsize%
  \chapter*{\listtablename}
   %\vskip 1\baselineskip
  \tdmpage
  {
    \setlength{\parskip}{0pt}%
      % Interligne simple pour la table des matières
      \renewcommand{\baselinestretch}{1}\small\normalsize%
    \@starttoc{lot}%
    }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numérotation des tableaux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand*\l@table[2]{%
  \settowidth{\largNum}{\hss #2}
  \ifnum \c@tocdepth >\m@ne%
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@%
    \setlength\@tempdima{3.4em}%ici pour changer l'espace entre le numero et le titre du tableau !! (défaut 2.2em)
        \noindent%
    \begingroup
      \pretolerance=10000
      \parindent \z@ \rightskip \tocRightMargin%
      \parfillskip -\tocRightMargin%
      \leavevmode \normalsize%
     \advance\leftskip \tocLeftMarginSecondLineTable %5.2em%8.4em%
      \hskip -\leftskip%
      {\tablename\mbox{\hspace{4pt}}#1}\nobreak%
       \leaders\hbox{$\m@th%
        \mkern \@dotsep mu\hbox{.}\mkern \@dotsep%
        mu$} \hfil\nobreak\hb@xt@%
                    \largNum{\hss #2}\par%
      \penalty\@highpenalty%
    \endgroup
  \fi}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liste des extraits de code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\lstlistoflistings{%
%    \renewcommand{\baselinestretch}{1}\small\normalsize%
    \chapter*{\lstlistlistingname}%
     %\vskip 1\baselineskip
    \tdmpage
    {
      \setlength{\parskip}{0pt}%
      \renewcommand{\baselinestretch}{1}\small\normalsize%
      \@starttoc{lol}%
    }
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numérotation des extraits de code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand*\l@lstlisting[2]{%
  \settowidth{\largNum}{\hss #2}
  \ifnum \c@tocdepth >\m@ne%
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@%
    \setlength\@tempdima{3.4em}%
        \noindent%
    \begingroup
      \pretolerance=10000
      \parindent \z@ \rightskip \tocRightMargin%
      \parfillskip -\tocRightMargin%
     \leavevmode \normalsize%
     \advance\leftskip \tocLeftMarginSecondLineFigure%
      \hskip -\leftskip%
      {\lstlistingname\mbox{\hspace{4pt}}#1}\nobreak%
       \leaders\hbox{$\m@th%
        \mkern \@dotsep mu\hbox{.}\mkern \@dotsep%
        mu$} \hfil\nobreak\hb@xt@%
                    \largNum{\hss #2}\par%
      \penalty\@highpenalty%
    \endgroup
  \fi}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liste des symboles basé sur l'environnement ``description''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newlength{\descril}
\setlength{\descril}{2cm}
\newenvironment{listofsymbols}[1][\descril]{%
  \cleardoublepage
  \chapter*{List of Symbols}
  \begin{list}{\parbox{}}{%
      \let\makelabel\losymblabel%
      \leftmargin #1%
      \labelwidth \leftmargin%
      \labelsep    1mm
      }}%
  {\end{list}}
\newcommand*\losymblabel[1]{\hspace\labelsep
  \normalfont\parbox{ #1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Liste des symboles basé sur l'environnement ``description''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newenvironment{listofabbr}[1][\descril]{%
 \cleardoublepage
  \chapter*{List of Abbreviations}
  \begin{list}{\mbox{}}{%
      \let\makelabel\losymblabel%
      \leftmargin #1%
      \labelwidth \leftmargin%
      \labelsep    1mm
      }}%
  {\end{list}}
\renewcommand*\losymblabel[1]{\hspace\labelsep
  \normalfont\mbox{ #1}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mise en forme des légendes pour les figures et tableaux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newlength{\largetiqfig}
\setlength{\largetiqfig}{2.2cm}
\newlength{\larglegendfig}
\setlength{\larglegendfig}{\textwidth-\largetiqfig}


%\long\def\@caption#1[#2]#3{\par\addcontentsline{\csname ext@#1\endcsname}{#1}
%{%
%    \protect\numberline{\csname the#1\endcsname}{\ignorespaces #2}}%
%  \begingroup\@parboxrestore\normalsize%
%  \pretolerance=10000
%  \ifx\figurestring#1
%     \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}{#1}\par%
%  \else
%     \ifx\tablestring#1
%        \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}{#1}\par%
%        \footnotesize \centering
%     \else
%        \@makecaption{\csname fnum@#1\endcsname}{\ignorespaces #3}\par%
%        \footnotesize \centering
%     \fi
%  \fi
%  \pretolerance=100
%  \endgroup}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\newlength{\tempabovecaptionskip}   %% Pour enlever l'éspace après le caption du tableau


%%%%
% On utilise le package caption pour définir le format des légendes
%
\captionsetup[figure]{font=singlespacing,labelfont=normalfont, justification=centering,labelsep=quad}
\captionsetup[table]{font=singlespacing,labelfont=normalfont, labelsep=quad,justification=centering,skip=1\baselineskip}
\captionsetup[algorithm]{font=singlespacing,labelfont=normalfont, labelsep=quad,justification=centering,skip=1\baselineskip}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Etant donne que la commande \captionsetup{list=no} n a pas marche pour enlever les legendes
% des figures et des tableaux des listes de figures et tableau, je triche en declarant un nouveau
% float pour les annexes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\DeclareCaptionType{figureap}[Figure-A]
\captionsetup[figureap]{font=singlespacing,justification=centering,labelsep=quad}
\renewenvironment{figureap}[1][]{%
  \vspace{\parskip}%   % ajout decalage vers bas
  \@float{figureap}%
  \centering
}{%
  \end@float
  \vskip -\parskip
  \vspace{-\parskip} % supprime decalage vers bas ajouté au dessus
}
\renewcommand{\thefigureap}{{\Roman{chapter}-\arabic{figureap}}}

\newsubfloat{figureap} %% Definition des subfigures en annexe


%\DeclareCaptionType{tableap}[Tableau-A]
%\captionsetup[tableap]{labelfont=normalfont,labelsep=quad,justification=centering,skip=2\baselineskip}
%\renewenvironment{tableap}[1][]{%
%  \vspace{\baselineskip}
%  \@float{tableap}%
%  \vspace{-\parskip}      
%  \vspace{-\baselineskip}
%  \centering
%}{%
%  \end@float
%  \vskip -\parskip % pareil que pour figure
%}
%\renewcommand{\thetableap}{{\Roman{chapter}-\arabic{tableap}}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numérotation des annexes en chiffres arabes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand{\appendix}{
  \par
%  \setcounter{chapter}{0}%
%  \setcounter{section}{0}%
%  \gdef\@chapapp{\appendixname}%
%  \gdef\thechapter{\@Roman\c@chapter}

  % On ajoute un nouveau compteur (incremente par makechapterhead) 
  % pour ne pas avoir a remettre chapter a 0 
  % car sans ca y probleme avec lien hyperref
  \setcounter{annexnum}{0}
%  \setcounter{chapter}{1}%
  \setcounter{chapter}{0}%
  \setcounter{section}{0}%
  \gdef\@chapapp{\appendixname}%
  \gdef\thechapter{\@Roman\c@annexnum}
  


  %
  % Empeche les sections et sous sections d'etre dans la table des matières
  % 

%  \addtocontents{toc}{\protect\setcounter{tocdepth}{0}}
%  \renewcommand{\thesection}{\arabic{section}}
  
  \newcommand{\numero}{\thesection}
  \newcommand{\changenumero}[1]{\renewcommand{\numero}{##1}}

  \renewcommand \thesection{\@arabic\c@section} % Enleve le numero du chapitre courant

  \def\sectionannexe{\@ifstar\unnumberedsection\numberedsection}

  \def\numberedsection{\@ifnextchar[%]
    \numberedsectionwithtwoarguments\numberedsectionwithoneargument}
  \def\unnumberedsection{\@ifnextchar[%]
    \unnumberedsectionwithtwoarguments\unnumberedsectionwithoneargument}
  \def\numberedsectionwithoneargument##1{\numberedsectionwithtwoarguments[##1]{##1}}
  \def\unnumberedsectionwithoneargument##1{\unnumberedsectionwithtwoarguments[##1]{##1}}
  \def\numberedsectionwithtwoarguments[##1]##2{%
    \ifhmode\par\fi
    \removelastskip
    \vskip .7ex \goodbreak
    \noindent
    \begingroup
      \leavevmode\Large\bfseries\raggedright
      \normalsize\bfseries
      \numero\quad
      ##2
      \par
    \endgroup
    \vskip .7ex \nobreak
  }
  \def\unnumberedsectionwithtwoarguments[##1]##2{%
    \ifhmode\par\fi
    \removelastskip
    \vskip .7ex\goodbreak
    \noindent
    \begingroup
      \leavevmode\normalfont\Large\bfseries\raggedright
      \leavevmode\normalfont\Large\bfseries\raggedright
      \normalsize\bfseries
      ##2
      \par
    \endgroup
    \vskip .7ex\nobreak
  }
  \def\section{\refstepcounter{section}\changenumero{\thesection} \sectionannexe}
  \def\subsection{\refstepcounter{subsection}\changenumero{\thesubsection} \sectionannexe}
  \def\subsubsection{\refstepcounter{subsubsection}\changenumero{\thesubsubsection} \sectionannexe}
  \def\paragraph{\refstepcounter{paragraph}\changenumero{\theparagraph} \sectionannexe}
  \def\subparagraph{\refstepcounter{subparagraph}\changenumero{\thesubparagraph} \sectionannexe}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Pour une bonne mise en page des annexes
  %% ANNEXE doit apparaitre dans la table des matieres lors de sa premiere apparition
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  \newcommand{\nameannexe}{\appendixname}
  \newcommand{\multiannexe}{
    \renewcommand{\nameannexe}{ANNEXES}
    
    \widowpenalty=10000
\clubpenalty=10000
\raggedbottom

  }

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Modification des annexes numérotées        % NOUVELLE DEMANDE : PAS REPETER LE TITRE SUR LA PAGE SUIVANTE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  \renewcommand{\@makechapterhead}[1]{{
        \normalfont\bfseries%
        \begin{center}
          \MakeUppercase{\appendixname~\thechapter} % ANNEXE 1
          \vskip 3\p@%
          \MakeUppercase{##1} % Title of annexe
        \end{center}
        \thispagestyle{empty}%

  }}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Modification des annexes sans numéro (ne devrait pas être utilisé)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  \renewcommand{\@makeschapterhead}[1]{{
    \thispagestyle{empty}
        \normalfont\bfseries%
        \begin{center}
           \MakeUppercase{##1} % TITLE OF ANNEXE
        \end{center}
        \vskip 5\p@%
        \thispagestyle{empty}
  }}

%  \renewenvironment{figure}{\begin{figureap}}{\end{figureap}}
%  \renewenvironment{table}{\begin{tableap}}{\end{tableap}}

}

%
% Position correcte de la numéroration des pages: 1' du côté, 0.75. du haut

\pagestyle{fancy}
\lhead[\thepage]{}
\chead{}
\rhead[]{\thepage}
\lfoot{}
\cfoot{}
\rfoot{}
\topmargin=  -\fanctmargin
\headsep=\fanchsep
\renewcommand{\headrulewidth}{0pt}

%
% Pour avoir des citations à interligne simple
%
\renewenvironment{quote}{
  \linespread{1}\small\normalsize\list{}{\rightmargin\leftmargin}%
  \item\relax
}{\endlist \linespread{2} \small \normalsize\noindent\hspace{-0.333em}\vskip 12\p@}

%
% Bonne présentation des figures et tables
%

\renewenvironment{figure}[1][]{%
  %~\\ % Saut de ligne avant figure 
  % espace ajouter directement  dans le texte (a l endroit ou se trouve la figure)
  \vspace{\parskip}%   % ajout decalage vers bas
  \@float{figure}%
  \centering
}{%
  \end@float
  % enleve le parskip ajouter a l endroit ou se trouve la figure dans le text (pas tjrs autour de la figure)
  \vskip -\parskip
  \vspace{-\parskip} % supprime decalage vers bas ajouté au dessus
}

\renewenvironment{table}[1][]{%
  %~\\ % Saut de ligne avant figure 
  \vspace{\baselineskip}
  % espace ajouter directement  dans le texte (a l endroit ou se trouve la figure)
  \@float{table}%
  % espace ajouter autour de la table (la ou se trouve la figure dans le document)
  \vspace{-\parskip}      
  \vspace{-\baselineskip}
  \centering
}{%
  \end@float
  \vskip -\parskip % pareil que pour figure
}

%
% Positionnement des flottants
%   on privilégie htbp (here, top, bottom, next page)
%
\def\fps@figure{htbp}
\def\fps@table{htbp}

%
% Énumération avec des lettres, aligné à gauche
% Modifi'e par LUIS (defs. de \par*) 
%
\renewcommand\enumerate{%
  \ifnum \@enumdepth >\thr@@\@toodeep\else
    \advance\@enumdepth\@ne
    \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \expandafter
      \list
        \csname label\@enumctr\endcsname
        {\usecounter
            {\@enumctr}% 
            \itemsep 0pt%
            \parsep \parskip% 
            \setlength{\partopsep}{0ex}%
           \def\makelabel##1{\hspace{-6pt}\hspace\labelsep ##1 \hspace\labelsep \hspace\labelsep}%
        }%
  \fi} 
\renewcommand{\labelenumi}{\alph{enumi}.}


% Pas de décallage dans l'itemize non plus, désactivation du custom itemize de babel french

%\frenchbsetup{StandardLists=true}

\let\tempitemize=\itemize
\renewcommand\itemize{%
 \setlength{\partopsep}{0pt}%
 \tempitemize
 \setlength{\itemsep}{0pt}%
 \def\makelabel##1{\hspace{-6pt}\hspace\labelsep ##1 \hspace\labelsep \hspace\labelsep} }


%
% Bibliographie
% Pas utilisé -Stéphane Mercier

\renewcommand{\bibname}{BIBLIOGRAPHY}

\theoremstyle{plain}
\newtheorem{definition}{Definition}[chapter]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rédéfinition de \cleardoublepage pour l'impression recto-verso
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\makeatletter
\def\cleardoublepage{\clearpage\if@twoside \ifodd\c@page\else
    \hbox{}
    \thispagestyle{empty}
    \newpage
    \if@twocolumn\hbox{}\newpage\fi\fi\fi}
\makeatother \clearpage{\pagestyle{empty}\cleardoublepage}