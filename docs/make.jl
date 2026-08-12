using Documenter, Pyrometers
mathengine = Documenter.MathJax3()
makedocs(
                sitename = "Pyrometers.jl",
                highlightsig = false,
                checkdocs = :none,
                format=Documenter.HTML(size_threshold = 2000 * 2^10 , 
                 mathengine = mathengine
                ),
                pages=[
                        #"Pyrometers"=>"index.md"
                        "Pyrometers" => "pyrometers.md",
                        "Examples"=>"pluto.md"
                        
                ]
               #
        )
#deploydocs(;
#         repo="github.com/Manarom/BandPyrometry"
#)