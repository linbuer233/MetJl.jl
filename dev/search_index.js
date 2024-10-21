var documenterSearchIndex = {"docs":
[{"location":"91-developer/#dev_docs","page":"Developer documentation","title":"Developer documentation","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"note: Contributing guidelines\nIf you haven't, please read the Contributing guidelines first.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If you want to make contributions to this package that involves code, then this guide is for you.","category":"page"},{"location":"91-developer/#First-time-clone","page":"Developer documentation","title":"First time clone","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: If you have writing rights\nIf you have writing rights, you don't have to fork. Instead, simply clone and skip ahead. Whenever upstream is mentioned, use origin instead.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If this is the first time you work with this repository, follow the instructions below to clone the repository.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fork this repo\nClone your repo (this will create a git remote called origin)\nAdd this repo as a remote:\ngit remote add upstream https://github.com/linbuer233/MetJl.jl","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"This will ensure that you have two remotes in your git: origin and upstream. You will create branches and push to origin, and you will fetch and update your local main branch from upstream.","category":"page"},{"location":"91-developer/#Linting-and-formatting","page":"Developer documentation","title":"Linting and formatting","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Install a plugin on your editor to use EditorConfig. This will ensure that your editor is configured with important formatting settings.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We use https://pre-commit.com to run the linters and formatters. In particular, the Julia code is formatted using JuliaFormatter.jl, so please install it globally first:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # Press ]\npkg> activate\npkg> add JuliaFormatter","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To install pre-commit, we recommend using pipx as follows:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"# Install pipx following the link\npipx install pre-commit","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"With pre-commit installed, activate it as a pre-commit hook:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit install","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To run the linting and formatting manually, enter the command below:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"pre-commit run -a","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Now, you can only commit if all the pre-commit tests pass.","category":"page"},{"location":"91-developer/#Testing","page":"Developer documentation","title":"Testing","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"As with most Julia packages, you can just open Julia in the repository folder, activate the environment, and run test:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"julia> # press ]\npkg> activate .\npkg> test","category":"page"},{"location":"91-developer/#Working-on-a-new-issue","page":"Developer documentation","title":"Working on a new issue","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"We try to keep a linear history in this repo, so it is important to keep your branches up-to-date.","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Fetch from the remote and fast-forward your local main\ngit fetch upstream\ngit switch main\ngit merge --ff-only upstream/main\nBranch from main to address the issue (see below for naming)\ngit switch -c 42-add-answer-universe\nPush the new local branch to your personal remote repository\ngit push -u origin 42-add-answer-universe\nCreate a pull request to merge your remote branch into the org main.","category":"page"},{"location":"91-developer/#Branch-naming","page":"Developer documentation","title":"Branch naming","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"If there is an associated issue, add the issue number.\nIf there is no associated issue, and the changes are small, add a prefix such as \"typo\", \"hotfix\", \"small-refactor\", according to the type of update.\nIf the changes are not small and there is no associated issue, then create the issue first, so we can properly discuss the changes.\nUse dash separated imperative wording related to the issue (e.g., 14-add-tests, 15-fix-model, 16-remove-obsolete-files).","category":"page"},{"location":"91-developer/#Commit-message","page":"Developer documentation","title":"Commit message","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Use imperative or present tense, for instance: Add feature or Fix bug.\nHave informative titles.\nWhen necessary, add a body with details.\nIf there are breaking changes, add the information to the commit message.","category":"page"},{"location":"91-developer/#Before-creating-a-pull-request","page":"Developer documentation","title":"Before creating a pull request","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"tip: Atomic git commits\nTry to create \"atomic git commits\" (recommended reading: The Utopic Git History).","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Make sure the tests pass.\nMake sure the pre-commit tests pass.\nFetch any main updates from upstream and rebase your branch, if necessary:\ngit fetch upstream\ngit rebase upstream/main BRANCH_NAME\nThen you can open a pull request and work with the reviewer to address any issues.","category":"page"},{"location":"91-developer/#Building-and-viewing-the-documentation-locally","page":"Developer documentation","title":"Building and viewing the documentation locally","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Following the latest suggestions, we recommend using LiveServer to build the documentation. Here is how you do it:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Run julia --project=docs to open Julia in the environment of the docs.\nIf this is the first time building the docs\nPress ] to enter pkg mode\nRun pkg> dev . to use the development version of your package\nPress backspace to leave pkg mode\nRun julia> using LiveServer\nRun julia> servedocs()","category":"page"},{"location":"91-developer/#Making-a-new-release","page":"Developer documentation","title":"Making a new release","text":"","category":"section"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"To create a new release, you can follow these simple steps:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Create a branch release-x.y.z\nUpdate version in Project.toml\nUpdate the CHANGELOG.md:\nRename the section \"Unreleased\" to \"[x.y.z] - yyyy-mm-dd\" (i.e., version under brackets, dash, and date in ISO format)\nAdd a new section on top of it named \"Unreleased\"\nAdd a new link in the bottom for version \"x.y.z\"\nChange the \"[unreleased]\" link to use the latest version - end of line, vx.y.z ... HEAD.\nCreate a commit \"Release vx.y.z\", push, create a PR, wait for it to pass, merge the PR.\nGo back to main screen and click on the latest commit (link: https://github.com/linbuer233/MetJl.jl/commit/main)\nAt the bottom, write @JuliaRegistrator register","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"After that, you only need to wait and verify:","category":"page"},{"location":"91-developer/","page":"Developer documentation","title":"Developer documentation","text":"Wait for the bot to comment (should take < 1m) with a link to a RP to the registry\nFollow the link and wait for a comment on the auto-merge\nThe comment should said all is well and auto-merge should occur shortly\nAfter the merge happens, TagBot will trigger and create a new GitHub tag. Check on https://github.com/linbuer233/MetJl.jl/releases\nAfter the release is create, a \"docs\" GitHub action will start for the tag.\nAfter it passes, a deploy action will run.\nAfter that runs, the stable docs should be updated. Check them and look for the version number.","category":"page"},{"location":"95-reference/#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"95-reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Pages = [\"95-reference.md\"]","category":"page"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [MetJl]","category":"page"},{"location":"95-reference/#MetJl.compute-Tuple{Any}","page":"Reference","title":"MetJl.compute","text":"hi = hello_world() A simple function to return \"Hello, World!\"\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.hello_world-Tuple{}","page":"Reference","title":"MetJl.hello_world","text":"hi = hello_world() A simple function to return \"Hello, World!\"\n\n\n\n\n\n","category":"method"},{"location":"95-reference/","page":"Reference","title":"Reference","text":"Modules = [calc]","category":"page"},{"location":"95-reference/#MetJl.calc.LCL-Tuple{Union{Unitful.Quantity{T, 𝐌 𝐋^-1 𝐓^-2, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐌 𝐋^-1 𝐓^-2, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝚯, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝚯, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝚯, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝚯, U}} where {L, S}} where {T, U}}","page":"Reference","title":"MetJl.calc.LCL","text":"LCL(pressure::Unitful.Pressure, T_air::Unitful.Temperature, dewpoint::Unitful.Temperature)\n\nCalculates the lifting condensation level (LCL) of a parcel of air.\n\nThe LCL is the pressure level at which a parcel of air becomes saturated and begins to condense when lifted dry adiabatically from the surface. This function uses an iterative approach to find the LCL by comparing the specific humidity of the parcel at a given pressure level to the specific humidity at the surface.\n\nArguments:\n\npressure: The pressure at the surface in Unitful pressure units.\nT_air: The temperature at the surface in Unitful temperature units.\ndewpoint: The dew point temperature at the surface in Unitful temperature units.\n\nReturns:\n\nA tuple of two Unitful quantities: the temperature and pressure at the LCL.\n\nExample:\n\njulia> calc.LCL(1013.25u\"hPa\", 25u\"°C\", 20u\"°C\")\n(291.9970210003215 K, 942 hPa)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.add_height_to_pressure-Tuple{Any, Any}","page":"Reference","title":"MetJl.calc.add_height_to_pressure","text":"addheightto_pressure(pressure, height)\n\nAdjusts the pressure level by adding a specified height, assuming standard atmospheric conditions.\n\nThis function uses a standard atmospheric model to convert the input pressure to a height, adds the specified height to this pressure level, and then converts the new height back to a pressure level using the same standard model.\n\nArguments:\n\npressure: The pressure level in a Unitful pressure unit.\nheight: The height to be added in a Unitful length unit.\n\nReturns:\n\nUnitful.Pressure: The new pressure level after adding the specified height, in hPa.\n\nExample:\n\njulia> calc.add_height_to_pressure(1000u\"hPa\",800u\"m\")\n908.303207892483 hPa\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.add_pressure_to_height-Tuple{Any, Any}","page":"Reference","title":"MetJl.calc.add_pressure_to_height","text":"addpressureto_height(height, pressure)\n\nAdjusts a given height by a certain pressure difference to find the new height at which the adjusted pressure would occur.\n\nThis function uses standard atmospheric pressure to height conversion functions to determine the new height. It first converts the input height to meters and pressure to hectopascals (hPa), then finds the pressure level corresponding to the input height. Finally, it subtracts the given pressure from this pressure level to find the new height at which this adjusted pressure would occur.\n\nArguments:\n\nheight: The initial height in a Unitful length unit.\npressure: The pressure difference to be subtracted in a Unitful pressure unit.\n\nReturns:\n\nA Unitful.Length: The new height at which the adjusted pressure would occur.\n\nExample:\n\njulia> add_pressure_to_height(1000u\"m\", 10u\"hPa\")\n1.0911907787454445 km\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.bulk_shear-Tuple{}","page":"Reference","title":"MetJl.calc.bulk_shear","text":"垂直风切\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.cape_cin-Tuple{Vector{Unitful.Quantity{Int64, 𝚯, Unitful.FreeUnits{(K,), 𝚯, nothing}}}, Vector{Unitful.Quantity{Int64, 𝚯, Unitful.FreeUnits{(K,), 𝚯, nothing}}}, Vector{Unitful.Quantity{Int64, 𝐌 𝐋^-1 𝐓^-2, Unitful.FreeUnits{(hPa,), 𝐌 𝐋^-1 𝐓^-2, nothing}}}, Vector{Unitful.Quantity{Int64, 𝚯, Unitful.FreeUnits{(K,), 𝚯, nothing}}}}","page":"Reference","title":"MetJl.calc.cape_cin","text":"capecin(Tair, dewpoint, pressure, T_stra)\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.dcape-Tuple{}","page":"Reference","title":"MetJl.calc.dcape","text":"downdraft_cape\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.dewpoint-Tuple{Any}","page":"Reference","title":"MetJl.calc.dewpoint","text":"dewpoint(vapor_pressure)\n\nCalculates the dew point temperature from the given vapor pressure.\n\nThis function uses the formula to convert the vapor pressure to the dew point temperature in Celsius, and then converts the result to Kelvin.\n\nArguments:\n\nvapor_pressure: The vapor pressure in a Unitful pressure unit.\n\nReturns:\n\nUnitful.Temperature: The dew point temperature in Kelvin.\n\nExample:\n\njulia> calc.dewpoint(6.11u\"hPa\")\n273.1455029235403 K\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.dewpoint_from_rh-Tuple{Union{Unitful.Quantity{T, 𝚯, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝚯, U}} where {L, S}} where {T, U}, Real}","page":"Reference","title":"MetJl.calc.dewpoint_from_rh","text":"dewpointfromrh(T_air::Unitful.Temperature, RH::Real)\n\nCalculates the dew point temperature given the air temperature and relative humidity.\n\nThis function uses the Magnus-Tetens approximation to estimate the saturation vapor pressure and iteratively adjusts the dew point temperature until the calculated vapor pressure matches the saturation vapor pressure corresponding to the dew point.\n\nArguments:\n\nT_air::Unitful.Temperature: The current air temperature.\nRH::Real: The relative humidity as a real number between 0 and 1.\n\nReturns:\n\nUnitful.Temperature: The dewpoint temperature in Kelvin.\n\nExample:\n\njulia> calc.dewpoint_from_rh(300u\"K\",0.5)\n288.701999999987 K\n\njulia> calc.dewpoint_from_rh([300,300,300]*u\"K\",[0.5,0.5,0.5])\n3-element Vector{Quantity{Float64, 𝚯, Unitful.FreeUnits{(K,), 𝚯, nothing}}}:\n 288.701999999987 K\n 288.701999999987 K\n 288.701999999987 K\n\nNote:\n\nThe function assumes that the input temperature is in a Unitful temperature unit and the relative humidity is a real number between 0 and 1.\nThe function uses an iterative approach to find the dew point temperature, which may not converge for very low temperatures or high humidities.\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.dewpoint_from_sh-Tuple{Any, Any}","page":"Reference","title":"MetJl.calc.dewpoint_from_sh","text":"dewpointfromsh(pressure, q)\n\nCalculates the dew point temperature from the given pressure and specific humidity.\n\nThis function uses the saturation vapor pressure formula to calculate the dew point temperature. It first converts the input pressure to hectopascals (hPa) and the specific humidity to grams per kilogram (g/kg). Then, it calculates the vapor pressure using the relationship between specific humidity, vapor pressure, and pressure. Finally, it uses the dewpoint function to find the dew point temperature corresponding to the calculated vapor pressure.\n\nArguments:\n\npressure::Unitful.Pressure: The atmospheric pressure.\nq::Unitful.Quantity: The specific humidity as a ratio (g/kg).\n\nReturns:\n\nUnitful.Temperature: The dew point temperature.\n\nExample:\n\njulia> calc.dewpoint_from_sh(1013.25u\"hPa\", 0.01u\"g/kg\")\n212.51310407992628 K\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.div-Tuple{}","page":"Reference","title":"MetJl.calc.div","text":"散度\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.parcel_profile-Tuple{Vector{Unitful.Quantity{Int64, 𝐌 𝐋^-1 𝐓^-2, Unitful.FreeUnits{(hPa,), 𝐌 𝐋^-1 𝐓^-2, nothing}}}, Union{Unitful.Quantity{T, 𝚯, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝚯, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝚯, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝚯, U}} where {L, S}} where {T, U}}","page":"Reference","title":"MetJl.calc.parcel_profile","text":"TODO 计算层结曲线\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.q_from_dewpoint-Tuple{Any, Any}","page":"Reference","title":"MetJl.calc.q_from_dewpoint","text":"q_from_dewpoint(pressure, dewpoint)\n\nCalculates the specific humidity (q) given the atmospheric pressure and dew point temperature.\n\nThis function uses the dew point temperature to calculate the saturation vapor pressure (e) using the vapor_pressure_from_dewpoint function, and then computes the specific humidity using the formula:\n\nq = 0.622 * e / (pressure - 0.378 * e)\n\nwhere:\n\nq is the specific humidity (in kg/kg),\ne is the saturation vapor pressure (in hPa),\npressure is the atmospheric pressure (in hPa),\ndewpoint is the dew point temperature (in Kelvin).\n\nArguments:\n\npressure::Unitful.Pressure: The atmospheric pressure.\ndewpoint::Unitful.Temperature: The dew point temperature.\n\nReturns:\n\nq Unitful.Quantity: The specific humidity (in g/kg).\n\nExample:\n\njulia> calc.q_from_dewpoint(983u\"hPa\", 298.15u\"K\")\n20.303761468426334 g kg⁻¹\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.q_from_vapor_pressure-Tuple{Any, Any}","page":"Reference","title":"MetJl.calc.q_from_vapor_pressure","text":"q_from_vapor_pressure(pressure, vapor_pressure)\n\nCalculates the specific humidity (q) from the given atmospheric pressure and vapor pressure.\n\nThis function uses the formula to convert vapor pressure to specific humidity, which is the ratio of the mass of water vapor to the total mass of the air. The formula takes into account the molecular weight difference between water vapor and dry air, and the assumption that the air behaves as an ideal gas.\n\nArguments:\n\npressure::Unitful.Pressure: The atmospheric pressure in which the vapor pressure is measured.\nvapor_pressure::Unitful.Pressure: The vapor pressure of water in the air.\n\nReturns:\n\nq Unitful.Quantity: The specific humidity in grams per kilogram (g/kg).\n\nExample: ```julia julia> calc.qfromvapor_pressure(1013.25u\"hPa\", 2.3u\"hPa\") 1.4131049133102709 g kg⁻¹\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.rh_from_dewpont-Tuple{Any, Any}","page":"Reference","title":"MetJl.calc.rh_from_dewpont","text":"rhfromdewpoint(T_air, dewpoint)\n\nCalculates the relative humidity (RH) given the air temperature and dew point temperature.\n\nThis function uses the formula for saturation vapor pressure as a function of temperature to calculate the vapor pressures at the air temperature and dew point temperature, then computes the relative humidity as the ratio of these two vapor pressures.\n\nArguments:\n\nT_air: The current air temperature in a Unitful temperature unit.\ndewpoint: The dew point temperature in a Unitful temperature unit.\n\nReturns:\n\nA Unitful.Quantity: The relative humidity as a real number between 0 and 1.\n\nExample:\n\njulia> calc.rh_from_dewpont(300u\"K\",288u\"K\")\n0.4779147464796312\n\njulia> calc.rh_from_dewpont([300,300,300]u\"K\",[288,288,288]u\"K\")\n3-element Vector{Float64}:\n 0.4779147464796312\n 0.4779147464796312\n 0.4779147464796312\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.theta-Tuple{Any, Any}","page":"Reference","title":"MetJl.calc.theta","text":"theta(T_air, Pressure)\n\nCalculates the potential temperature (θ) given the air temperature and pressure.\n\nThis function uses the formula for potential temperature, which is a measure of the temperature of air after it has been brought adiabatically to a reference pressure, usually 1000 hPa. The formula is:\n\nθ = T_air * (1000 / Pressure)^0.286\n\nArguments:\n\nT_air: The current air temperature in a Unitful temperature unit.\nPressure: The atmospheric pressure in a Unitful pressure unit.\n\nReturns:\n\nA Unitful.Quantity: The potential temperature in Kelvin.\n\nExample:\n\njulia> calc.theta(25u\"°C\", 800u\"hPa\")\n317.79793975452117 K\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.vapor_pressure_from_dewpoint-Tuple{Any}","page":"Reference","title":"MetJl.calc.vapor_pressure_from_dewpoint","text":"vaporpressurefrom_dewpoint(dewPoint)\n\nCalculates the saturation vapor pressure from the dew point temperature.\n\nThis function uses the Tetens formula to estimate the saturation vapor pressure from the dew point temperature. The Tetens formula is an approximation that is valid for temperatures between 0°C and 50°C.\n\nArguments:\n\ndewPoint: The dew point temperature in a Unitful temperature unit.\n\nReturns:\n\nUnitful.Pressure: The saturation vapor pressure in hectopascals (hPa).\n\nExample:\n\njulia> calc.vapor_pressure_from_dewpoint(20u\"°C\")\n23.397012933719008 hPa\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.vo-Tuple{}","page":"Reference","title":"MetJl.calc.vo","text":"涡度\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.w_speed-Tuple{}","page":"Reference","title":"MetJl.calc.w_speed","text":"垂直速度\n\n\n\n\n\n","category":"method"},{"location":"95-reference/#MetJl.calc.wind_chill_index-Tuple{Any, Any}","page":"Reference","title":"MetJl.calc.wind_chill_index","text":"windchillindex(windspeed, Tair)\n\nCalculates the wind chill index (WCI) based on the wind speed and air temperature.\n\nThis function implements the formula for the wind chill index, which estimates the perceived temperature felt by the human body due to the cooling effect of the wind. The formula is:\n\nWCI = (10 * sqrt(windspeed) - windspeed + 10.5) * (33 - T_air)\n\nArguments:\n\nwind_speed: The wind speed in a Unitful velocity unit.\nT_air: The air temperature in a Unitful temperature unit.\n\nReturns:\n\nA Unitful.Quantity: The wind chill index with units of kg*cal/m^2/h.\n\nExample:\n\njulia> calc.wind_chill_index(20u\"m/s\", 290u\"K\")\n568.824956732432 cal kg m⁻² s⁻¹\n\n\n\n\n\n","category":"method"},{"location":"90-contributing/#contributing","page":"Contributing guidelines","title":"Contributing guidelines","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"First of all, thanks for the interest!","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"We welcome all kinds of contribution, including, but not limited to code, documentation, examples, configuration, issue creating, etc.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"Be polite and respectful, and follow the code of conduct.","category":"page"},{"location":"90-contributing/#Bug-reports-and-discussions","page":"Contributing guidelines","title":"Bug reports and discussions","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you think you found a bug, feel free to open an issue. Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.","category":"page"},{"location":"90-contributing/#Working-on-an-issue","page":"Contributing guidelines","title":"Working on an issue","text":"","category":"section"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If you found an issue that interests you, comment on that issue what your plans are. If the solution to the issue is clear, you can immediately create a pull request (see below). Otherwise, say what your proposed solution is and wait for a discussion around it.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"tip: Tip\nFeel free to ping us after a few days if there are no responses.","category":"page"},{"location":"90-contributing/","page":"Contributing guidelines","title":"Contributing guidelines","text":"If your solution involves code (or something that requires running the package locally), check the developer documentation. Otherwise, you can use the GitHub interface directly to create your pull request.","category":"page"},{"location":"","page":"MetJl","title":"MetJl","text":"CurrentModule = MetJl","category":"page"},{"location":"#MetJl","page":"MetJl","title":"MetJl","text":"","category":"section"},{"location":"","page":"MetJl","title":"MetJl","text":"Documentation for MetJl.","category":"page"},{"location":"#Contributors","page":"MetJl","title":"Contributors","text":"","category":"section"},{"location":"","page":"MetJl","title":"MetJl","text":"<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->\n<!-- prettier-ignore-start -->\n<!-- markdownlint-disable -->\n<table>\n    <tbody>\n        <tr>\n            <td align=\"center\" valign=\"top\" width=\"14.28%\">\n                <a href=\"http://linbuer233.github.io\">\n                    <img src=\"https://avatars.githubusercontent.com/u/94788509?v=4?s=100\" width=\"100px;\"\n                        alt=\"Linbuer233\" />\n                    <br />\n                    <sub>\n                        <b>Linbuer233</b>\n                    </sub>\n                </a>\n                <br />\n                <a href=\"#code-abelsiqueira\" title=\"Code\">💻</a>\n                <a href=\"#projectManagement-abelsiqueira\" title=\"Project Management\">📆</a>\n                <a href=\"#doc-abelsiqueira\" title=\"Documentation\">📖</a>\n                <a href=\"#maintenance-abelsiqueira\" title=\"Maintenance\">🚧</a>\n                <a href=\"#review-tmigot\" title=\"Reviewed Pull Requests\">👀</a>\n            </td>\n            <td align=\"center\" valign=\"top\" width=\"14.28%\">\n                <a href=\"https://github.com/wzzcharles\">\n                    <img src=\"https://avatars.githubusercontent.com/u/102020995?v=4?s=100\" width=\"100px;\"\n                        alt=\"wzzcharles\" />\n                    <br />\n                    <sub>\n                        <b>wzzcharles</b>\n                    </sub>\n                </a>\n                <br />\n                <a href=\"#code-tmigot\" title=\"Code\">💻</a>\n                <a href=\"#doc-tmigot\" title=\"Documentation\">📖</a>\n                <a href=\"#review-tmigot\" title=\"Reviewed Pull Requests\">👀</a>\n                <a href=\"#ideas-oxinabox\" title=\"Ideas, Planning, & Feedback\">🤔</a>\n            </td>\n        </tr>\n    </tbody>\n</table>\n<!-- markdownlint-restore -->\n<!-- prettier-ignore-end -->\n\n<!-- ALL-CONTRIBUTORS-LIST:END -->","category":"page"},{"location":"#Special-Thanks","page":"MetJl","title":"Special Thanks","text":"","category":"section"},{"location":"","page":"MetJl","title":"MetJl","text":"Our project stands on the shoulders of giants. We are grateful to the open source community for the tools, libraries, and frameworks that have made our work easier and more efficient.","category":"page"},{"location":"","page":"MetJl","title":"MetJl","text":"<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->\n<!-- prettier-ignore-start -->\n<!-- markdownlint-disable -->\n<table>\n    <tbody>\n        <tr>\n            <td align=\"center\" valign=\"top\" width=\"14.28%\">\n                <a href=\"https://abelsiqueira.com\">\n                    <img src=\"https://avatars.githubusercontent.com/u/1068752?v=4?s=100\" width=\"100px;\"\n                        alt=\"Abel Soares Siqueira\" />\n                    <br />\n                    <sub><b>Abel Soares Siqueira</b>\n                    </sub>\n                </a>\n            </td>\n        </tr>\n    </tbody>\n</table>\n<!-- markdownlint-restore -->\n<!-- prettier-ignore-end -->\n\n<!-- ALL-CONTRIBUTORS-LIST:END -->","category":"page"}]
}
