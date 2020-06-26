\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{phd}[20/02/2018 v1.0]

\LoadClass[twoside,fontsize=24pt,a4paper,openright]{book}	

%\usepackage{jsfds}

%\RequirePackage[utf8]{inputenc}	
\usepackage{scrextend} %Forcer la 4eme  de couverture en page pair
\usepackage{minitoc}
\RequirePackage[round, sort]{natbib}
\RequirePackage[english, frenchb]{babel}
\usepackage[fixlanguage]{babelbib}
\RequirePackage{ifpdf, color}
\RequirePackage{textcomp}
\usepackage[utf8]{inputenc}
\RequirePackage[OT1]{fontenc}
%\RequirePackage{lmodern}
\RequirePackage{tikz}
\RequirePackage{fix-cm}
\RequirePackage{graphicx}
\RequirePackage{sectsty}
\RequirePackage{shorttoc}
\RequirePackage{fancyhdr}
\RequirePackage{pifont}
\RequirePackage{stmaryrd}
\RequirePackage{bm}
\RequirePackage{caption}
\RequirePackage{breakcites}
\RequirePackage{amsmath,amsthm}
\RequirePackage{helvet}
\RequirePackage{eurosym}
\RequirePackage{microtype}
\RequirePackage{amssymb}
\RequirePackage{fancybox}
\RequirePackage{blindtext}
%\RequirePackage[]{algorithm2e}
\RequirePackage[]{algorithm}
\RequirePackage{rotating}
\RequirePackage{dsfont}
\RequirePackage{frcursive}
\RequirePackage{appendix}
\RequirePackage{enumerate, enumitem}
\renewcommand{\normalsize}{\fontsize{11}{12}\selectfont}

\setlist[description]{%
  topsep=10pt,               
  itemsep=10pt,               
  font={\normalfont}, 
}

%\RequirePackage{algorithm}
\RequirePackage{algorithmic}

\usepackage{xstring}
\usepackage{changepage}
\RequirePackage{multicol}
\strictpagecheck

%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Define your number of chapters here!!! %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand\ChapNum{4}

% Graphiques et couleurs pour jsfds
\RequirePackage{xcolor}
\RequirePackage{graphicx}
% Compteurs de chapitres


%%%%%% Debut code onglets
\RequirePackage{background}
\usetikzlibrary{positioning,calc}
%

\newcommand\DeactivateBG{\backgroundsetup{contents={}}}
\newcommand\ActivateBG{\backgroundsetup{
scale=1,
angle=0,
opacity=0.5,
contents={
\ifnum\value{page}>0%
\ifnum\value{page}<121%
\checkoddpage
\ifoddpage
\begin{tikzpicture}[remember picture,overlay]
\node[rotate=90,minimum width=\paperheight,minimum height=0.5cm,fill=gray!40,font=\large\bfseries,text=white,anchor=south west] at ($(current page.south east)$) {};
\node[rotate=90,minimum width=1/\ChapNum*\paperheight,minimum height=0.5cm,fill=\BoxColor,font=\large\bfseries,text=white,anchor=south west] at ([yshift=(\ChapNum-(\thechapter))/\ChapNum*\paperheight]$(current page.south east)$) 
{\ifthenelse{\equal{\thechapter}{4}}{Conclusion}{Chapter \thechapter}};
%{\ifthenelse{\equal{\thechapter}{1}}{Introduction}{\ifthenelse{\equal{\thechapter}{4}}{Conclusion}{Chapter \thechapter}}};
\end{tikzpicture}
\else
\begin{tikzpicture}[remember picture,overlay]
\node[rotate=90,minimum width=\paperheight,minimum height=0.5cm,fill=gray!40,font=\large\bfseries,text=white,anchor=south west] at ([xshift=-\paperwidth+0.5cm]$(current page.south east)$) {};
\node[rotate=90,minimum width=1/\ChapNum*\paperheight,minimum height=0.5cm,fill=\BoxColor,font=\large\bfseries,text=white,anchor=south west] at ([yshift=\ChapNum-(\thechapter))/\ChapNum*\paperheight,xshift=-\paperwidth+0.5cm]$(current page.south east)$) {\ifthenelse{\equal{\thechapter}{4}}{Conclusion}{Chapter \thechapter}};%{\ifthenelse{\equal{\thechapter}{1}}{Introduction}{\ifthenelse{\equal{\thechapter}{4}}{Conclusion}{Chapter \thechapter}}};
\end{tikzpicture}
\fi
\fi
\fi
}}}

\definecolor{bleuvert}{RGB}{102, 255, 204}
\definecolor{Prune}{RGB}{99,0,60}
\definecolor{indianred2}{RGB}{238,99,99}
\definecolor{orange2}{RGB}{238,154,0}
\definecolor{steelblue4}{RGB}{54,100,139}
\definecolor{turquoise3}{RGB}{0,197,205}
\definecolor{steelblue4}{RGB}{54,100,139}
\definecolor{mint}{RGB}{76, 137, 37}
 \definecolor{greenAgro}{HTML}{028887}
\definecolor{tryanpurple}{RGB}{136,2,70}
\newcommand\BoxColor{%
  \ifcase\thechapter  greenAgro!100!black\or  mint!100!black\or orange2!85!black\or steelblue4!100!black\or Prune!100!black\or black!40!black\or red!40!black\or blue!40!black\else green!40!black\fi}

\sectionfont{\color{\BoxColor}\sffamily\bfseries}
\subsectionfont{\sffamily \mdseries}
\subsubsectionfont{\sffamily\bfseries}

%%%% End onglets
\definecolor{rosefonce}{RGB}{197, 88, 88}
\definecolor{bleuvertgtis}{RGB}{88, 88, 197}



\RequirePackage{atveryend}
\RequirePackage{titletoc}
 \usepackage{tocloft}
%\renewcommand\cftchapafterpnum{\vskip6pt}
%\renewcommand\cftsecafterpnum{\vskip6pt}
%\renewcommand\cftsubsecafterpnum{\vskip6pt}
%
%% Adjust sectional unit title fonts in ToC
%\renewcommand{\cftchapfont}{\sffamily}
%\renewcommand{\cftsecfont}{\sffamily}
%\renewcommand{\cftsubsecfont}{\sffamily}
\renewcommand{\cfttoctitlefont}{\sffamily\Huge\bfseries}
\contentsmargin[1cm]{0cm}
\titlecontents{chapter}[0em]{\vskip12pt\bfseries\sffamily}
{\thecontentslabel\enspace}
{\hspace{1.05em}}
{ \hfill\contentspage}[\vskip 6pt]

\titlecontents{section}[1em]{\sffamily}
{\thecontentslabel\enspace}
{}
{\titlerule*[1pc]{.}\quad\contentspage}[\vskip 0pt]

\titlecontents{subsection}[2.7em]{\sffamily}
{\thecontentslabel\enspace}
{}
{\titlerule*[1pc]{.}\quad\contentspage}[\vskip 3pt]

%\usepackage{etoolbox}
%\pretocmd{\contentsname}{\sffamily}{}{}
%\RequirePackage[inner=30mm,outer=20mm,top=30mm,bottom=30mm]{geometry}
\RequirePackage[inner=40mm,outer=30mm,top=40mm,bottom=40mm]{geometry}

\RequirePackage{setspace} 

\RequirePackage{hyperref} 
%\hypersetup{colorlinks=true, citecolor=blue, filecolor=black, linkcolor=black, urlcolor=black}
\hypersetup{
    bookmarks=true,         % show bookmarks bar?
    unicode=false,          % non-Latin characters in Acrobat’s bookmarks
    pdftoolbar=true,        % show Acrobat’s toolbar?
    pdfmenubar=true,        % show Acrobat’s menu?
    pdffitwindow=false,     % window fit to page when opened
    pdfstartview={FitH},    % fits the width of the page to the window
    pdftitle={Network inference from incomplete abundance data},    % title
    pdfauthor={Raphaelle Momal},     % author
    pdfsubject={Applied mathematics},   % subject of the document
    pdfcreator={Raphaelle Momal},   % creator of the document
    pdfnewwindow=true,      % links in new PDF window
    colorlinks=true,       % false: boxed links; true: colored links
    linkcolor=black,%bleuvertgtis,          % color of internal links (change box color with linkbordercolor)
    citecolor=blue2,        % color of links to bibliography
    filecolor=rosefonce,      % color of file links
    urlcolor=blue1           % color of external links
}
 
\definecolor{blue1}{RGB}{55,126,184}
\definecolor{blue2}{RGB}{29,120,166}%49,113,165
\RequirePackage{framed} 
\RequirePackage{color} 
\usetikzlibrary{calc,shapes,backgrounds,arrows,automata,shadows,positioning}
\newcommand{\EqRef}[1] {Equation (\ref{#1})}
\onehalfspacing
\RequirePackage{ifthen} 
\RequirePackage{csquotes}
%\frenchspacing
%\FrenchFootnotes

%%%%%%%%%%%%%% En-tête de chapitre %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\makeatletter
%\def\thickhrulefill{\leavevmode \leaders \hrule height 1ex \hfill \kern \z@}
%\def\@makechapterhead#1{%
%	\vspace*{-30\p@}%
%	{\parindent \z@ \raggedleft \reset@font
%				\scshape {\color{\BoxColor}{\@chapapp{} {\fontsize{100pt}{12pt}\selectfont{\textbf{\thechapter}}}}}
%			\par\nobreak
%			\interlinepenalty\@M
%		\Huge \bfseries {\color{\BoxColor}{#1}}\par\nobreak
%		%\vspace*{1\p@}%
%		\hrulefill
%		\par\nobreak
%		\vskip 10\p@
%	}}
%\def\@makeschapterhead#1{%
%	\vspace*{-10\p@}%
%		{\parindent \z@ \raggedleft \reset@font
%			\scshape \vphantom{\@chapapp{} \thechapter}
%			\par\nobreak
%			\interlinepenalty\@M
%		\Huge \bfseries #1 \par\nobreak
%		%\vspace*{1\p@}%
%		\hrulefill
%		\par\nobreak
%		\vskip 10\p@
%	}}
%
%
\makeatletter
\def\thickhrulefill{\leavevmode \leaders \hrule height 1ex \hfill \kern \z@}
\def\@makechapterhead#1{%
	\vspace*{-10\p@}%
	{
	\parindent \z@ \raggedright \reset@font
				\bfseries
				\ifthenelse{\equal{\thechapter}{0}}{}
				{\color{\BoxColor}{ {\fontsize{80pt}{12pt}\selectfont{\textsf{\textbf{\thechapter}}}}}}%\@chapapp{}
			\par\nobreak
			\interlinepenalty\@M
		\Huge  \sffamily \scshape{\color{\BoxColor}{#1}}\par\nobreak
		\vspace*{-10\p@}%
		{\color{\BoxColor}\hrulefill}
		\par\nobreak
		\vskip 10\p@
	}}
%\def\@makeschapterhead#1{%
%	\vspace*{-10\p@}%
%		{\parindent \z@ \raggedleft \reset@font
%			\scshape \vphantom{\@chapapp{} \thechapter}
%			\par\nobreak
%			\interlinepenalty\@M
%		\Huge \bfseries #1 \par\nobreak
%		%\vspace*{1\p@}%
%		\hrulefill
%		\par\nobreak
%		\vskip 10\p@
%	}}

%%%%%%%%%%%%%%%%%%%% Fin en-tête de chapitre %%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Tête de parties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\makeatletter
\newcommand*{\jolipart}[1]{%
	\begin{center}
		\begin{Huge}\color{black}
			#1
		\end{Huge}
		\end{center}
		\vspace{1cm}
	\begin{center}
%\begin{Huge}
		%\ding{167}
	%\end{Huge}
	\hrulefill
	\end{center}
	}


\parttitlefont{\jolipart}

\partnumberfont{\sc\normalsize}



%%%%%%%%%%%%%%%% COULEURS grises %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\definecolor{grisfonce}{RGB}{49,49,49}
\definecolor{grisclair}{RGB}{111,111,111}
\definecolor{blanc}{RGB}{255,255,255}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% RÉGLAGE DE LA BIBLIOGRAPHIE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ATTENTION, on utilise Biblatex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%\RequirePackage[style=verbose-trad1,hyperref,backend=biber]{biblatex} %Pour le style 
%% biblatex avec gestion des mots clés des entrées bibliographiques
%% + les abréviations idem, ibidem, op cit. Nécessite de prévoir un fichier
%%  bib contenant ces mots-clés. Vous pouvez changer l'option verbose-trad1 
%% pour une autre (cf. manuel de biblatex).
%% ATTENTION : l'option backend=biber est utilisée. En effet les dernières
%% versions de biblatex préconisent l'option biber pour profiter des dernières
%% nouveautés. Si vous choisissez cette option (backend=biber), il vous faut
%% simplement compiler avec $biber nomdufichier (sans le .tex). Si vous voulez
%% en rester à bibtex,il vous suffit de placer backend=bibtex.
%
%\bibliography{biblio} %Nom du fichier bibtex a utiliser.
%\defbibheading{primary}{\section*{Sources primaires}} % Catégorie (sous forme de section) pour le tri de la bibliographie.
%
%\defbibheading{secondary}{\section*{Sources secondaires}} % Catégorie (sous forme de section) pour le tri de la bibliographie.
%
%\renewcommand*{\newunitpunct}{\addcomma\space} % Virgule entre les parties d'une référence (merci a Josquin Debaz)
%
%%%%% Ajustements pour la bibliographie
%%\DeclareFieldFormat[article]{volume}{\textbf{#1}}  %Le numéro de volume en gras
%\DeclareFieldFormat[article]{number}{\textit{#1}} %Le numéro dans le volume en italique
%%\DeclareFieldFormat{pages}{page(s): #1} % page(s) en toutes lettres, si on veut...
%
%\renewcommand{\mkibid}[1]{\emph{#1}} %Les locutions latines en italique (comme ibid, loc.cit. , etc.) Merci à Josquin Debaz
%
%% Et pour mettre le in en italique dans la référence biblio (merci encore à Josquin Debaz)
%\DefineBibliographyStrings{french}{%
%	in = {\emph{in}}%
%}


%%%%%%%%%%%%%%%%%% Interligne simple pour citations et vers %%%%%%%%%
%% Environnements: QUOTE, QUOTATION, VERSE
% Copyright (C) 2002 Mike Nolta <mrnolta@princeton.edu>, GPL V. 2.0 and
% later version. http://www.physics.princeton.edu/cosmology/computing/PrincetonThesis.cls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\let\orig@quote\quote
\let\endorig@quote\endquote
\renewenvironment*{quote}
	{\begin{orig@quote}\begin{singlespace}}
	{\end{singlespace}\end{orig@quote}}

\let\old@quotation\quotation
\let\endold@quotation\endquotation
\renewenvironment*{quotation}
	{\begin{old@quotation}\begin{singlespace}}
		{\end{singlespace}\end{old@quotation}}

\let\old@verse\verse
\let\endold@verse\endverse
\renewenvironment*{verse}
	{\begin{old@verse}\begin{singlespace}}
		{\end{singlespace}\end{old@verse}}

%%%%%%%%%%%%%%%% TITRE THESE PAGE DE GARDE %%%%%%%%%%%%%%%%

\newcommand*{\logoEntreprise}[1]{\gdef\@logoEntreprise{#1}}
\logoEntreprise{}

\newcommand*{\NNTnumber}[1]{\gdef\@NNTnumber{#1}}
\NNTnumber{}

\newcommand*{\scalelogoEntreprise}[1]{\gdef\@scalelogoEntreprise{#1}}
\scalelogoEntreprise{}

%logo université

\newcommand*{\logouniversite}[1]{\gdef\@logouniversite{#1}}
\logouniversite{}

\newcommand*{\scalelogouniversite}[1]{\gdef\@scalelogouniversite{#1}}
\scalelogouniversite{}

% Logo Labo
\newcommand*{\logolabo}[1]{\gdef\@logolabo{#1}}
\logolabo{}

\newcommand*{\scalelogolabo}[1]{\gdef\@scalelogolabo{#1}}
\scalelogolabo{}

% Nom de l'université
\newcommand*{\nomdeuniversite}[1]{\gdef\@nomdeuniversite{#1}}
\nomdeuniversite{}

%%%% Éventuel sous-titre
\newcommand{\lesoustitre}[1]{\gdef\@lesoustitre{#1}}
\lesoustitre{}

%%%% Discipline
\newcommand{\discipline}[1]{\gdef\@discipline{#1}}
\discipline{}

%%%% Nom du directeur de thèse
\newcommand{\dirdethese}[1]{\gdef\@dirdethese{#1}}
\dirdethese{}

%%%% Titre du directeur de thèse
\newcommand{\titredudirdethese}[1]{\gdef\@titredudirdethese{#1}}
\titredudirdethese{}

%%%% Jury
\newcommand{\jury}[1]{\gdef\@jury{#1}}
\jury{}

%%%% Laboratoire / départment / unité
\newcommand{\unite}[1]{\gdef\@unite{#1}}
\unite{}

%%%%% École doctorale
\newcommand{\ecoledoc}[1]{\gdef\@ecoledoc{#1}}
\ecoledoc{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mise en page de la page de garde
\makeatletter
\def\maketitle{%
	\thispagestyle{empty}
	\clearpage
	\begin{tikzpicture}[remember picture,overlay,line width=0mm]
	\draw [draw=white,fill=white]
		(current page.north west) rectangle (\paperwidth,1);
	\node[xshift=0.085\paperwidth,yshift=1.5cm,text=white,font=\bf\Large] {
	\includegraphics[scale=\@scalelogouniversite]{\@logouniversite}
	};
	\node[xshift=.6\paperwidth,yshift=1.5cm,text=white,font=\bf\Large] {
	\includegraphics[scale=\@scalelogolabo]{\@logolabo}
	};
	\node[yshift=-0.5cm] at (current page.center)
{\parbox{15cm}{
	\begin{center}
		\sffamily\Large Thèse présentée pour obtenir le  grade de docteur \\ \@nomdeuniversite
	\end{center}
\begin{center}
\begin{huge}
\color{greenAgro}\ding{91}
\end{huge}
\end{center}
\begin{center}
\@unite \\
\@ecoledoc
\end{center}
	\begin{center}
		\textbf{Discipline : \@discipline}
	\end{center}
\vspace{0.5cm}
\noindent\rule{15cm}{0.5mm}
	\begin{center}
		{\setlength{\baselineskip}{1.4\baselineskip}
		\begin{Large}
			\textbf{\@title}
		\end{Large}
		\par}
	\end{center}
%	\begin{center}
%		
%		\textbf{\@lesoustitre}
%	\end{center}
\noindent\rule{15cm}{0.5mm}
	\begin{flushright}
		\textsc{Par : \textbf{\@author}}
	\end{flushright}
\vspace{0.5cm}
%\bigskip
	\begin{flushleft}
		\textbf{Sous la direction de} \textsc{\@dirdethese}, \@titredudirdethese
	\end{flushleft}
	\textsc{\underline{Membres du jury:}}
	\vspace{0.5cm}
	\begin{normalsize}
		\@jury
	\end{normalsize}
	\vspace{1.5cm}
	\begin{flushright}
		\textbf{Date de soutenance :} \@date
	\end{flushright}
}
};
\node[xshift=0cm,yshift=-23.5cm,text=white,font=\bf\Large] {
	\includegraphics[scale=\@scalelogoEntreprise]{\@logoEntreprise}
	};
\node[xshift=-3cm,yshift=1.2cm] at ($(current page.south east)$) {
	\begin{footnotesize}
		\sffamily NNT=\@NNTnumber
	\end{footnotesize}
	};
\end{tikzpicture}
}

\makeatother


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% En-têtes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% Pour supprimer les en-têtes et pied de page gênants par exemple juste avant un chapitre sur une page de droite
\newcommand{\clearemptydoublepage}{%
\newpage{\pagestyle{empty}\cleardoublepage}}
%%%% ... et utiliser la commande \clearemptydoublepage juste avant le \chapter


\fancyhf{}	% on annule le fancy automatique


%%%%%%%%%% Gérer les en-têtes dans les frontmatter mainmatter et backmatter

\RequirePackage{etoolbox}

\appto\frontmatter{\pagestyle{fancy}
\renewcommand{\sectionmark}[1]{}
\renewcommand{\chaptermark}[1]{\markboth{\textit{#1}}{}}
\fancyhead[LE,RO]{\small\thepage}
	\fancyhead[RO]{\small\leftmark}		% \rightmark = section courante
	\fancyhead[LE]{\small\leftmark}		% \leftmark = chapitre courant
	\fancyfoot[C]{\thepage}				% marque la page au centre
}

\appto\mainmatter{\pagestyle{plain}%\pagestyle{fancy}
%\renewcommand{\sectionmark}[1]{\markright{\textit{\thesection.\ #1}}}
%\renewcommand{\chaptermark}[1]{\markboth{\chaptername~\thechapter~--\ \textit{#1}}{}}
%\fancyhead[LE,RO]{\small\thepage}
%	\fancyhead[RO]{\small\rightmark}					% \rightmark = section courante
%	\fancyhead[LE]{\small \leftmark}	% \leftmark = chapitre courant
%	\fancyfoot[C]{\thepage}								% marque la page au centre
}
%\appto\mainmatter{\pagestyle{fancy}
%\renewcommand{\headrulewidth}{0pt}
%\fancyhead{}
%	\fancyfoot[C]{\thepage}				% marque la page au centre
%}
\appto\backmatter{\pagestyle{fancy}
\renewcommand{\sectionmark}[1]{\markright{\thesection.\ #1}}
\renewcommand{\chaptermark}[1]{\markboth{\chaptername~\thechapter~--\ #1}{}}
\fancyhead[LE,RO]{\small\thepage}
	\fancyhead[RO]{}		% \rightmark = section courante
	\fancyhead[LE]{}		% \leftmark = chapitre courant
	\fancyfoot[C]{\thepage}	% marque la page au centre
}

%%%%%%%%%%%%% Appendices après chaque chapitre %%%%%%%%%%%%%%%%
% Start of subappendices environment
\AtBeginEnvironment{subappendices}{%
\chapter*{Appendices}
\addcontentsline{toc}{chapter}{Appendices}
\counterwithin{figure}{section}
\counterwithin{table}{section}
}

% End of subappendices environment
\AtEndEnvironment{subappendices}{%
\counterwithout{figure}{section}
\counterwithout{table}{section}
}
%%%%%%%%%%%%%% FAIRE DES ENCARTS %%%%%%%%%%%%%%%%%%%%%%%%%%
% Couleurs :

\definecolor{fondtitre}{RGB}{85,85,85}
\definecolor{fonddeboite}{RGB}{232,232,232}
\definecolor{shadecolor}{RGB}{232,232,232}

\ProcessOptions

%%%%%%%%%% Boitemagique  %%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand*{\boitemagique}[2]{
\begin{center}
\begin{tikzpicture}
% la boite
\node[rectangle,draw=fondtitre!100,fill=gray!40,inner sep=10pt, inner ysep=20pt] (mabox)
{
\begin{minipage}{12cm}
#2
\end{minipage}
};
% le titre de la boite
\node[fill=\BoxColor, text=white, rectangle] at (mabox.north) {\sffamily\textbf{#1}};
\end{tikzpicture}
\end{center}
}



%%%%%%%%%% petiteboite  %%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand*{\petiteboite}[2]{
\begin{center}
\begin{tikzpicture}
% la boite
\node[rectangle,draw=fondtitre!100,fill=gray!15,inner sep=10pt, inner ysep=20pt] (mabox)
{
\begin{minipage}{13cm}
#2
\end{minipage}
};
% le titre de la boite
\node[fill=\BoxColor!60, text=\BoxColor!20, rectangle] at (mabox.north) {\sffamily\textbf{#1}};
\end{tikzpicture}
\end{center}
}

%%%%%% boitesimple %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand*{\boitesimple}[1]{%
\begin{center}
\begin{minipage}{12cm}
\begin{shaded}
	#1
\end{shaded}
\end{minipage}
\end{center}
}

%%%%%%%%% New command for the thesis %%%%%%%%%%%


\newcommand{\pkg}[1]{{\fontseries{b}\selectfont #1}}
\newcommand{\proglang}[1]{{\fontseries{b}\selectfont #1}}
\providecommand{\keywords}[1]{\textbf{\textit{Keywords---}} #1}
\newcommand{\code}[1]{\texttt{#1}}

\newcommand{\WhitePage}
{
	\newpage
	\clearemptydoublepage
	\strut
	\newpage
}

%%%%%% Pour le JSS

\usepackage[]{graphicx}\usepackage[]{color}
\makeatletter
\def\maxwidth{ %
  \ifdim\Gin@nat@width>\linewidth
    \linewidth
  \else
    \Gin@nat@width
  \fi
}
\makeatother

\definecolor{fgcolor}{rgb}{0.345, 0.345, 0.345}
\newcommand{\hlnum}[1]{\textcolor[rgb]{0.686,0.059,0.569}{#1}}%
\newcommand{\hlstr}[1]{\textcolor[rgb]{0.192,0.494,0.8}{#1}}%
\newcommand{\hlcom}[1]{\textcolor[rgb]{0.678,0.584,0.686}{\textit{#1}}}%
\newcommand{\hlopt}[1]{\textcolor[rgb]{0,0,0}{#1}}%
\newcommand{\hlstd}[1]{\textcolor[rgb]{0.345,0.345,0.345}{#1}}%
\newcommand{\hlkwa}[1]{\textcolor[rgb]{0.161,0.373,0.58}{\textbf{#1}}}%
\newcommand{\hlkwb}[1]{\textcolor[rgb]{0.69,0.353,0.396}{#1}}%
\newcommand{\hlkwc}[1]{\textcolor[rgb]{0.333,0.667,0.333}{#1}}%
\newcommand{\hlkwd}[1]{\textcolor[rgb]{0.737,0.353,0.396}{\textbf{#1}}}%
\let\hlipl\hlkwb

\usepackage{framed}
\makeatletter
\newenvironment{kframe}{%
 \def\at@end@of@kframe{}%
 \ifinner\ifhmode%
  \def\at@end@of@kframe{\end{minipage}}%
  \begin{minipage}{\columnwidth}%
 \fi\fi%
 \def\FrameCommand##1{\hskip\@totalleftmargin \hskip-\fboxsep
 \colorbox{shadecolor}{##1}\hskip-\fboxsep
     % There is no \\@totalrightmargin, so:
     \hskip-\linewidth \hskip-\@totalleftmargin \hskip\columnwidth}%
 \MakeFramed {\advance\hsize-\width
   \@totalleftmargin\z@ \linewidth\hsize
   \@setminipage}}%
 {\par\unskip\endMakeFramed%
 \at@end@of@kframe}
\makeatother

\definecolor{shadecolor}{rgb}{.97, .97, .97}
\definecolor{messagecolor}{rgb}{0, 0, 0}
\definecolor{warningcolor}{rgb}{1, 0, 1}
\definecolor{errorcolor}{rgb}{1, 0, 0}
\newenvironment{knitrout}{}{} % an empty environment to be redefined in TeX

\usepackage{alltt}

\newcommand\abstractname{\textcolor{green}{Abstract}}
\if@titlepage
  \newenvironment{abstract}{%
      \titlepage
      \null\vfil
      \@beginparpenalty\@lowpenalty
      %\begin{center}%
        %\bfseries \abstractname
        \large \bfseries {\color{\BoxColor}{\abstractname}}
        \@endparpenalty\@M
      %\end{center}}%
     {\par\vfil\null\endtitlepage}
\else
  \newenvironment{abstract}{%
      \if@twocolumn
        \section*{\abstractname}%
      \else
        \small
       % \begin{center}% \bfseries {\color{\BoxColor}{#1}}
         % {\bfseries \abstractname\vspace{-.5em}}%
         {\large \bfseries {\color{\BoxColor}{\abstractname}}\vspace{-.5em}}
        %\end{center}%
        \quotation
      \fi}
      {\if@twocolumn\else\endquotation\fi}
\fi

%%

\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}

