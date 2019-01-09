[[AbstractFFTs]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "8d59c3b1463b5e0ad05a3698167f85fac90e184d"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "0.3.2"

[[ArnoldiMethod]]
deps = ["DelimitedFiles", "LinearAlgebra", "Random", "SparseArrays", "StaticArrays", "Test"]
git-tree-sha1 = "2b6845cea546604fb4dca4e31414a6a59d39ddcd"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.0.4"

[[Arpack]]
deps = ["BinaryProvider", "Libdl", "LinearAlgebra", "Random", "SparseArrays", "Test"]
git-tree-sha1 = "1ce1ce9984683f0b6a587d5bdbc688ecb480096f"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.3.0"

[[AxisAlgorithms]]
deps = ["Compat", "WoodburyMatrices"]
git-tree-sha1 = "99dabbe853e4f641ab21a676131f2cf9fb29937e"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "0.3.0"

[[AxisArrays]]
deps = ["Compat", "Dates", "IntervalSets", "IterTools", "Random", "RangeArrays", "Test"]
git-tree-sha1 = "2e2536e9e6f27c4f8d09d8442b61a7ae0b910c28"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.3.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BinDeps]]
deps = ["Compat", "Libdl", "SHA", "URIParser"]
git-tree-sha1 = "12093ca6cdd0ee547c39b1870e0c9c3f154d9ca9"
uuid = "9e28174c-4ba2-5203-b857-d8d62c4213ee"
version = "0.8.10"

[[BinaryProvider]]
deps = ["Libdl", "Pkg", "SHA", "Test"]
git-tree-sha1 = "055eb2690182ebc31087859c3dd8598371d3ef9e"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.3"

[[BufferedStreams]]
deps = ["Compat", "Test"]
git-tree-sha1 = "5d55b9486590fdda5905c275bb21ce1f0754020f"
uuid = "e1450e63-4bb3-523b-b2a4-4ffa8c0fd77d"
version = "1.0.0"

[[Cairo]]
deps = ["BinDeps", "Colors", "Compat", "Graphics", "Homebrew", "Libdl", "WinRPM"]
git-tree-sha1 = "a427098d5aa2808504d94b8ed9fc5740ceaf71d0"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "0.5.6"

[[CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays", "Test"]
git-tree-sha1 = "254cf73ea369d2e39bfd6c5eb27a2296cfaed68c"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.0"

[[CodecZlib]]
deps = ["BinaryProvider", "Libdl", "Test", "TranscodingStreams"]
git-tree-sha1 = "e3df104c84dfc108f0ca203fd7f5bbdc98641ae9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.5.1"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random", "Test"]
git-tree-sha1 = "f73b0e10f2a5756de7019818a41654686da06b09"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.7.5"

[[ColorVectorSpace]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "a890f08e61b40e9843d7177206da61229a3603c8"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.6.2"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "InteractiveUtils", "Printf", "Reexport", "Test"]
git-tree-sha1 = "9f0a0210450acb91c730b730a994f8eef1d3d543"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.9.5"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "ec61a16eed883ad0cfa002d7489b3ce6d039bb9a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "1.4.0"

[[Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Test", "UUIDs"]
git-tree-sha1 = "4f396f3aa55d3077a619daa058029cc297573bdd"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.7.2"

[[ComputationalResources]]
deps = ["Test"]
git-tree-sha1 = "89e7e7ed20af73d9f78877d2b8d1194e7b6ff13d"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.0"

[[Conda]]
deps = ["Compat", "JSON", "VersionParsing"]
git-tree-sha1 = "fb86fe40cb5b35990e368709bfdc1b46dbb99dac"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.1.1"

[[CoordinateTransformations]]
deps = ["Compat", "Rotations", "StaticArrays"]
git-tree-sha1 = "47f05d0b7f4999609f92e657147df000818c1f24"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.5.0"

[[CustomUnitRanges]]
deps = ["Test"]
git-tree-sha1 = "0a106457a1831555857e18ac9617279c22fc393b"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "0.2.0"

[[DataStructures]]
deps = ["InteractiveUtils", "OrderedCollections", "Random", "Serialization", "Test"]
git-tree-sha1 = "ca971f03e146cf144a9e2f2ce59674f5bf0e8038"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.15.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distances]]
deps = ["LinearAlgebra", "Printf", "Random", "Statistics", "Test"]
git-tree-sha1 = "0e37d8a95bafbeb9c800ef27ab6f443d29e17610"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.7.4"

[[Distributed]]
deps = ["LinearAlgebra", "Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["Distributed", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "c24e9b6500c037673f0241a2783472b8c3d080c7"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.16.4"

[[FFTViews]]
deps = ["CustomUnitRanges", "FFTW", "Test"]
git-tree-sha1 = "9d7993227ca7c0fdb6b31deef193adbba11c8f4e"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.2.0"

[[FFTW]]
deps = ["AbstractFFTs", "BinaryProvider", "Compat", "Conda", "Libdl", "LinearAlgebra", "Reexport", "Test"]
git-tree-sha1 = "29cda58afbf62f35b1a094882ad6c745a47b2eaa"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "0.2.4"

[[FileIO]]
deps = ["Pkg", "Random", "Test"]
git-tree-sha1 = "1a114d08094e7267ba8d4d684f930d5a722184de"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.0.4"

[[FixedPointNumbers]]
deps = ["Test"]
git-tree-sha1 = "b8045033701c3b10bf2324d7203404be7aef88ba"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.5.3"

[[Glob]]
deps = ["Compat", "Test"]
git-tree-sha1 = "c72f1fcb7d17426de1e8af2e948dfb3de1116eed"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.2.0"

[[Graphics]]
deps = ["Colors", "Compat", "NaNMath"]
git-tree-sha1 = "e3ead4211073d4117a0d2ef7d1efc5c8092c8412"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "0.4.0"

[[Gtk]]
deps = ["BinDeps", "Cairo", "Compat", "Graphics", "Homebrew", "Libdl", "Reexport", "Serialization", "Test", "WinRPM"]
git-tree-sha1 = "261bfcf4910db53990f289c1338aa8eb3b89fb44"
uuid = "4c0ca9eb-093a-5379-98c5-f87ac0bbbf44"
version = "0.16.4"

[[GtkReactive]]
deps = ["Cairo", "Colors", "Dates", "FixedPointNumbers", "Graphics", "Gtk", "IntervalSets", "Reactive", "Reexport", "RoundingIntegers", "Test"]
git-tree-sha1 = "6d5596de3fd5d6d2f012f6f36d3b226eac5e2b47"
uuid = "27996c0f-39cd-5cc1-a27a-05f136f946b6"
version = "0.5.3"

[[HTTPClient]]
deps = ["Compat", "LibCURL"]
git-tree-sha1 = "161d5776ae8e585ac0b8c20fb81f17ab755b3671"
uuid = "0862f596-cf2d-50af-8ef4-f2be67dfa83f"
version = "0.2.1"

[[Homebrew]]
deps = ["BinDeps", "InteractiveUtils", "JSON", "Libdl", "Test", "Unicode"]
git-tree-sha1 = "f01fb2f34675f9839d55ba7238bab63ebd2e531e"
uuid = "d9be37ee-ecc9-5288-90f1-b9ca67657a75"
version = "0.7.1"

[[IdentityRanges]]
deps = ["OffsetArrays", "Test"]
git-tree-sha1 = "f5ca23a08397288924a7535cd898d8ccbcde5ba5"
uuid = "bbac6d45-d8f3-5730-bfe4-7a449cd117ca"
version = "0.2.0"

[[ImageAxes]]
deps = ["AxisArrays", "Colors", "FixedPointNumbers", "ImageCore", "MappedArrays", "Reexport", "SimpleTraits", "Test"]
git-tree-sha1 = "5735ec90843acaa67a4624611921c686cdf4efbf"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.5.0"

[[ImageCore]]
deps = ["ColorTypes", "Colors", "FFTW", "FixedPointNumbers", "Graphics", "MappedArrays", "OffsetArrays", "PaddedViews", "Random", "Statistics", "Test"]
git-tree-sha1 = "5e7b1f49c80541860e08a7ea91805a24c1641f19"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.7.3"

[[ImageDistances]]
deps = ["Colors", "Distances", "LinearAlgebra", "ProgressMeter", "Test"]
git-tree-sha1 = "a5de7b61f6fa98fb93c39857fa43cf40ca383b28"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.1.1"

[[ImageFiltering]]
deps = ["CatIndices", "ColorVectorSpace", "Colors", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "FixedPointNumbers", "ImageCore", "LinearAlgebra", "Logging", "MappedArrays", "OffsetArrays", "Random", "StaticArrays", "Statistics", "Test", "TiledIteration"]
git-tree-sha1 = "9bb70bff423fb26b45817da5f127d71e662bfe50"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.5.2"

[[ImageMetadata]]
deps = ["AxisArrays", "ColorVectorSpace", "Colors", "Dates", "FixedPointNumbers", "ImageAxes", "ImageCore", "IndirectArrays", "Statistics", "Test"]
git-tree-sha1 = "e8763955937114ea15cf3bfb5ae00588e19b58c3"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.5.2"

[[ImageMorphology]]
deps = ["ImageCore", "Test"]
git-tree-sha1 = "e94f43b9ff76f3a3810bfdd9b3d2fbcacbc26fd0"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.1.1"

[[ImageShow]]
deps = ["Base64", "ColorTypes", "Colors", "FileIO", "FixedPointNumbers", "ImageCore", "OffsetArrays", "Requires", "Test"]
git-tree-sha1 = "98eb96a852fd2d6f0905cbe4d215ec2113805b46"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.1.2"

[[ImageTransformations]]
deps = ["AxisAlgorithms", "ColorTypes", "ColorVectorSpace", "Colors", "CoordinateTransformations", "FixedPointNumbers", "IdentityRanges", "ImageCore", "Interpolations", "LinearAlgebra", "OffsetArrays", "StaticArrays", "Test"]
git-tree-sha1 = "18ae1c0a8df31549a9452ceac93751d4aa166071"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.7.1"

[[ImageView]]
deps = ["AxisArrays", "Cairo", "ColorVectorSpace", "Colors", "FileIO", "FixedPointNumbers", "Graphics", "Gtk", "GtkReactive", "Images", "MappedArrays", "RoundingIntegers", "StatsBase", "Test"]
git-tree-sha1 = "be3e41d030ad6792225814844dbe77e085006f42"
uuid = "86fae568-95e7-573e-a6b2-d8a6b900c9ef"
version = "0.8.2"

[[Images]]
deps = ["AxisArrays", "Base64", "ColorTypes", "ColorVectorSpace", "Colors", "FileIO", "FixedPointNumbers", "Graphics", "ImageAxes", "ImageCore", "ImageDistances", "ImageFiltering", "ImageMetadata", "ImageMorphology", "ImageShow", "ImageTransformations", "IndirectArrays", "LinearAlgebra", "MappedArrays", "OffsetArrays", "Random", "Reexport", "SIUnits", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "Test", "TiledIteration"]
git-tree-sha1 = "01291372b2c3a67191a8cc2ff291f33f5485ac8b"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.17.0"

[[IndirectArrays]]
deps = ["Compat", "Test"]
git-tree-sha1 = "b6e249be10a3381b2c72ac82f2d13d70067cb2bd"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "0.5.0"

[[Inflate]]
deps = ["Pkg", "Printf", "Random", "Test"]
git-tree-sha1 = "b7ec91c153cf8bff9aff58b39497925d133ef7fd"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.1"

[[InteractiveUtils]]
deps = ["LinearAlgebra", "Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "Test", "WoodburyMatrices"]
git-tree-sha1 = "1ec5b78dbfb5ae2b3d7a4987a9287fed4cc31ded"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.11.1"

[[IntervalSets]]
deps = ["Compat"]
git-tree-sha1 = "9dc556002f23740de13946e8c2e41798e09a9249"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.3.1"

[[IterTools]]
deps = ["SparseArrays", "Test"]
git-tree-sha1 = "79246285c43602384e6f1943b3554042a3712056"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.1.1"

[[JLD2]]
deps = ["CodecZlib", "DataStructures", "FileIO", "LinearAlgebra", "Mmap", "Printf", "Random", "Test"]
git-tree-sha1 = "3ba90ff93e1d5b9b2103588051c2d349fae54dac"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.1.2"

[[JSON]]
deps = ["Dates", "Distributed", "Mmap", "Sockets", "Test", "Unicode"]
git-tree-sha1 = "1f7a25b53ec67f5e9422f1f551ee216503f4a0fa"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.20.0"

[[LibCURL]]
deps = ["BinaryProvider", "Compat", "Libdl", "Printf"]
git-tree-sha1 = "6339c87cb76923a3cf947fcd213cbc364355c9c9"
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.4.1"

[[LibExpat]]
deps = ["Compat"]
git-tree-sha1 = "fde352ec13479e2f90e57939da2440fb78c5e388"
uuid = "522f3ed2-3f36-55e3-b6df-e94fee9b0c07"
version = "0.5.0"

[[LibGit2]]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libz]]
deps = ["BufferedStreams", "Random", "Test"]
git-tree-sha1 = "d405194ffc0293c3519d4f7251ce51baac9cc871"
uuid = "2ec943e9-cfe8-584d-b93d-64dcb6d567b7"
version = "1.0.0"

[[LightGraphs]]
deps = ["ArnoldiMethod", "Base64", "DataStructures", "DelimitedFiles", "Distributed", "Inflate", "LinearAlgebra", "Markdown", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics", "Test"]
git-tree-sha1 = "c7222c370d5cf6d4e08ae40bddd8c0db6852dfb1"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.2.0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Compat"]
git-tree-sha1 = "c443e1c8d58a4e9f61b708ad0a88286c7042145b"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.4.4"

[[MappedArrays]]
deps = ["Test"]
git-tree-sha1 = "923441c5ac942b60bd3a842d5377d96646bcbf46"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.2.1"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Measures]]
deps = ["Test"]
git-tree-sha1 = "ddfd6d13e330beacdde2c80de27c1c671945e7d9"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.0"

[[MetaGraphs]]
deps = ["Base64", "JLD2", "LightGraphs", "Random", "Test"]
git-tree-sha1 = "82533eb68e57a1f39e4d371b43de24fd81a3b567"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.6.0"

[[Missings]]
deps = ["Dates", "InteractiveUtils", "SparseArrays", "Test"]
git-tree-sha1 = "adc26d2ee85a49c413464110d922cf21efc9d233"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.3.1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[NaNMath]]
deps = ["Compat"]
git-tree-sha1 = "ce3b85e484a5d4c71dd5316215069311135fa9f2"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.2"

[[OffsetArrays]]
deps = ["DelimitedFiles", "Test"]
git-tree-sha1 = "7d1442cb06fbfbc4fea936c3c56b38daffd22d3b"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "0.9.1"

[[OrderedCollections]]
deps = ["Random", "Serialization", "Test"]
git-tree-sha1 = "85619a3f3e17bb4761fe1b1fd47f0e979f964d5b"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.0.2"

[[PDMats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "b6c91fc0ab970c0563cbbe69af18d741a49ce551"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.9.6"

[[POMDPs]]
deps = ["Distributions", "LibGit2", "Pkg", "Random", "Statistics", "Test"]
git-tree-sha1 = "4214e0a48a4c746ecde390eed7bc0e2205dc40a0"
uuid = "a93abf59-7444-517b-a68a-c42f96afdd7d"
version = "0.7.2"

[[PaddedViews]]
deps = ["OffsetArrays", "Test"]
git-tree-sha1 = "7da3e7e1a58cffbf10177553ae95f17b92516912"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.4.2"

[[Parameters]]
deps = ["Markdown", "OrderedCollections", "REPL", "Test"]
git-tree-sha1 = "70bdbfb2bceabb15345c0b54be4544813b3444e4"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.10.3"

[[Pkg]]
deps = ["Dates", "LibGit2", "Markdown", "Printf", "REPL", "Random", "SHA", "UUIDs"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[ProgressMeter]]
deps = ["Distributed", "Printf", "Random", "Test"]
git-tree-sha1 = "48058bc11607676e5bbc0b974af79106c6200787"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "0.9.0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra", "Test"]
git-tree-sha1 = "3ce467a8e76c6030d4c3786e7d3a73442017cdc0"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.0.3"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RangeArrays]]
deps = ["Compat"]
git-tree-sha1 = "d925adfd5b01cb46fde89dc9548d167b3b136f4a"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.1"

[[Ratios]]
deps = ["Compat"]
git-tree-sha1 = "fd159bead0a24e6270fd0573a340312bd4645cc2"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.3.0"

[[Reactive]]
deps = ["DataStructures", "Distributed", "Test"]
git-tree-sha1 = "5862d915387ebb954016f50a88e34f79a9e5fcd2"
uuid = "a223df75-4e93-5b7c-acf9-bdd599c0f4de"
version = "0.8.3"

[[Reel]]
deps = ["Test", "VideoIO"]
git-tree-sha1 = "a6cab7f04d936324144228babf9aa0ca9a0d780a"
uuid = "71555da5-176e-5e73-a222-aebc6c6e4f2f"
version = "1.1.0"

[[Reexport]]
deps = ["Pkg"]
git-tree-sha1 = "7b1d07f411bc8ddb7977ec7f377b97b158514fe0"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "0.2.0"

[[Requires]]
deps = ["Test"]
git-tree-sha1 = "f6fbf4ba64d295e146e49e021207993b6b48c7d1"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "0.5.2"

[[Rmath]]
deps = ["BinaryProvider", "Libdl", "Random", "Statistics", "Test"]
git-tree-sha1 = "9a6c758cdf73036c3239b0afbea790def1dabff9"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.5.0"

[[Rotations]]
deps = ["LinearAlgebra", "Random", "StaticArrays", "Statistics", "Test"]
git-tree-sha1 = "40f63689787ad81cd014361405997b78a5595944"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "0.9.2"

[[RoundingIntegers]]
deps = ["Test"]
git-tree-sha1 = "293ba0ab32218b9ffd596040224228def84f8da0"
uuid = "d5f540fe-1c90-5db3-b776-2e2f362d9394"
version = "0.2.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SIUnits]]
deps = ["Compat", "TexExtensions"]
git-tree-sha1 = "224d83b62711fe7e429454aace2c97eb2cf36031"
uuid = "b9d75638-96e3-5676-bdf0-e9c958f63a55"
version = "0.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools", "Test"]
git-tree-sha1 = "c0a542b8d5e369b179ccd296b2ca987f6da5da0a"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.8.0"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures", "Random", "Test"]
git-tree-sha1 = "03f5898c9959f8115e30bc7226ada7d0df554ddd"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "0.3.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["BinDeps", "BinaryProvider", "Libdl", "Test"]
git-tree-sha1 = "0b45dc2e45ed77f445617b99ff2adf0f5b0f23ea"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "0.7.2"

[[StaticArrays]]
deps = ["InteractiveUtils", "LinearAlgebra", "Random", "Statistics", "Test"]
git-tree-sha1 = "1eb114d6e23a817cd3e99abc3226190876d7c898"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "0.10.2"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsBase]]
deps = ["DataStructures", "DelimitedFiles", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "Test"]
git-tree-sha1 = "7b596062316c7d846b67bf625d5963a832528598"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.27.0"

[[StatsFuns]]
deps = ["Rmath", "SpecialFunctions", "Test"]
git-tree-sha1 = "d14bb7b03defd2deaa5675646f6783089e0556f0"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.7.0"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Test]]
deps = ["Distributed", "InteractiveUtils", "Logging", "Random"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TexExtensions]]
deps = ["Compat"]
git-tree-sha1 = "092ad55ed044aa5ab31ee800d4ae5bec526a8f09"
uuid = "9b435220-3ad3-5d4f-b1ea-1e7b29ae9b13"
version = "0.1.0"

[[TiledIteration]]
deps = ["OffsetArrays", "Test"]
git-tree-sha1 = "58f6f07d3b54a363ec283a8f5fc9fb4ecebde656"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.2.3"

[[TranscodingStreams]]
deps = ["Pkg", "Random", "Test"]
git-tree-sha1 = "a34a2d588e2d2825602bf14a24216d5c8b0921ec"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.8.1"

[[URIParser]]
deps = ["Test", "Unicode"]
git-tree-sha1 = "6ddf8244220dfda2f17539fa8c9de20d6c575b69"
uuid = "30578b45-9adc-5946-b283-645ec420af67"
version = "0.4.0"

[[UUIDs]]
deps = ["Random"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Vec]]
deps = ["LinearAlgebra", "Printf", "Random", "StaticArrays"]
git-tree-sha1 = "e31c2b197a0a794e43e55e1b441ef1dc21d37444"
repo-rev = "master"
repo-url = "https://github.com/sisl/Vec.jl.git"
uuid = "44eeaf0b-fee4-471f-9310-ed6585cb3142"
version = "0.0.0"

[[VersionParsing]]
deps = ["Compat"]
git-tree-sha1 = "c9d5aa108588b978bd859554660c8a5c4f2f7669"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.1.3"

[[VideoIO]]
deps = ["BinDeps", "ColorTypes", "FileIO", "FixedPointNumbers", "Glob", "Homebrew", "ImageCore", "ImageView", "Libdl", "Requires", "Test"]
git-tree-sha1 = "6ec4f79ac04606587b0cad836c71acd51a30b169"
uuid = "d6d074c3-1acf-5d4c-9a43-ef38773959a2"
version = "0.4.0"

[[WinRPM]]
deps = ["BinDeps", "Compat", "HTTPClient", "LibExpat", "Libdl", "Libz", "URIParser"]
git-tree-sha1 = "2a889d320f3b77d17c037f295859fe570133cfbf"
uuid = "c17dfb99-b4f7-5aad-8812-456da1ad7187"
version = "0.4.2"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Test"]
git-tree-sha1 = "21772c33b447757ec7d3e61fcdfb9ea5c47eedcf"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.4.1"
