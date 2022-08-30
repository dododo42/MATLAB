classdef Mastermind < handle
    
%==========================================================================
% Bug:             Engine after move doesnt save the attempt (to fix)
%
% Version:         1.0.1 - BETA
%
% Description:     Run this class to play Mastermind on MATLAB. Drag the
%                  pegs in the position you want and base your next try on 
%                  the game feedback:
%                   - white peg: right color and wrong position
%                   - black peg: right color and right position
%                  Use the hint button to exploit an engine based on the
%                  Donald Knuth algorithm, which ensures a solution in 
%                  maximum 5 moves.
%               
% Author:          Eduardo Maria Polli
%                  eduardo.polli94@gmail.com
%               
% Date:            August 26, 2022
%
% Aknowledgments:  Thanks to Brian Moore for the awsome library on
%                  FileExchange!
%--------------------------------------------------------------------------

    properties (Access = private , Constant = true)
        
        % Board colors:
        windowColor = 0.75*ones(1,3);                                      % Background color                                
        lightColor = 0.9*ones(1,3);                                        % Light detail color
        darkColor = 0.45*ones(1,3);                                        % Dark detail color
                
        % Board dimensions:
        nRows = 10;                                                        % Number of attempts
        nCols = 4;                                                         % Code length
        resultsShift = 12;                                                 % Distance of result pegs [px]
        initialGap = 10;                                                   % Space left from border [px]
        pegRowGap = 12;                                                    % Vertical space between pegs [px]
        pegColGap = 24;                                                    % Horizontal space between pegs [px]
        borderWidth = 12;                                                  % Border width [px]
        detailLineWidth = 2;                                               % Width of detail boundary lines
                    
        % Peg appearance:
        pegColors = {'#ff0000','#00ff00','#0000ff',...                     % Available colors
                    '#ff7f00','#00ffff','#ffff00'};
        pegRadius = 25;                                                    % Peg radius [px]
        resultsRadius = 10;                                                % Feedback peg radius [px]
        lightCoeffH = .3;                                                  % Horizontal offset of peg light direction
        lightCoeffV = -.2;                                                 % Vertical offset of peg light direction
        
        % End Game Screen:
        resultBaseGap = 30;                                                % Text-box distance
        resultLowerGap = 120;                                              % Lower text-box distance
        
        % Engine:
        feedbackList = [0,0,0,0,0,1,1,1,1,2,2,2,3,4;...                    % All possible peg feedbecks for engine hint
                        0,1,2,3,4,0,1,2,3,0,1,2,0,0]';
        
    end
    
    properties (Access = private)
        
        % GUI:
        gameAxis;                                                          % Axis object
        button1Position;                                                   % Exit button handle
        button2Position;                                                   % Hint button handle
        button3Position;                                                   % Play button handle
        isData;                                                            % Bool for using data.mat splashStart arts
        
        % Game parameters.
        isPlaying;                                                         % Bool for ongoing game
        solution;                                                          % Numeric code to be unlocked
        currentRow;                                                        % Current attempt
        colorIdx;                                                          % Selected peg color                                                   
        defaultPegPosition;                                                % Saves initial peg position while dragging
        attempt;                                                           % Numeric code tried
        multipleGame;                                                      % Bool for multiple games in one sessions
        
        % Peg Coordinates:
        pegX;                                                              % Colored peg X coordinates
        pegY;                                                              % Colored peg Y coordinates
        pegHoleX;                                                          % Peg holes X coordinates
        pegHoleY;                                                          % Peg holes Y coordinates
        resultsX;                                                          % Result peg X coordinates
        resultsY;                                                          % Result peg Y coordinates
        
        % Peg graphics:
        pegGRX = struct;                                                   % Home screen splashStart art CData
        alphamap;                                                          % Home screen splashStart art AlphaMap
        
        % Peg objects:
        coloredPeg = struct;                                               % Colroed pegs handle
        holePeg = struct;                                                  % Peg holes handle
        resultsPeg = struct;                                               % Result pegs handle
        
        % splashStart arts:
        splashStart;                                                       % Start screen splashStart art handle
        splashEnd = struct;                                                         % End screen splashStart art handle
        GRXdata;                                                           % End screen splashStart art CData
        msgPosition;                                                       % End screen banner position
        
        % Engine:
        allSolutions;                                                      % List of all non-attempted solutions
        possibleSolutions;                                                 % List of possible solutions based on feedbacks
        solutionCount;                                                     % Count per color of the solution code
        
    end
    
    methods (Access = private)
        
        % =================================================================
        % GAME FUNCTIONS
        % =================================================================
        
        %% Start New Game:
        function newGame(app)
            
            % Prepare Board if this is the second game:
            if app.isData
                delete(app.splashStart)
                if app.multipleGame
                    delete(app.splashEnd.box)
                    delete(app.splashEnd.text)
                    for j = 1 : 5
                        if j < 5
                            delete(app.splashEnd.peg(j).obj)
                        end
                        delete(app.splashEnd.details(j).line)
                    end
                end
            end
            plotPegHoles(app);
            plotResultsPegs(app);
            plotMainPegs(app);
            
            % Reset game parameters:
            app.isPlaying = 1;
            app.currentRow = 1;
            app.attempt = 7*ones(1,4);
            app.solution = ceil(length(app.pegColors)*rand(1,app.nCols));
            
            % Activate mouse interactivity:
            set(gcf,'WindowButtonDownFcn',@app.determineAction)
            
            % Initialise Engine variables:
            app.solutionCount = countSolutionNumber(app,app.solution);
            app.allSolutions = getAllCombinations(app);
            app.possibleSolutions = app.allSolutions;
        end
        
        %% Number of pegs for each color in the solution code:
        function count = countSolutionNumber(app,x)
            
            % Count number of pegs in the solution code for each color
            count = zeros(1,length(app.pegColors));
            for j = 1 : length(app.pegColors)
                count(j) = sum(x==j);
            end
        end
        
        %% Decide the mouse callback based on the click coordinates:
        function determineAction(app,~,~)
            
            % Get click coordinates:
            mouseCoord = get(gca,'CurrentPoint');
            
            % Check cordinates position:
            if mouseCoord(1,2) < app.button3Position(2)+app.button3Position(4)
                if mouseCoord(1,1) > app.button1Position(1)
                    % Use buttons:
                    makeAction(app);
                elseif mouseCoord(1,1) < app.button1Position(1)
                    if app.isPlaying
                        % Move pegs:
                        determinePeg(app);
                    end
                end
            end
        end
        
        %% Determine which colored peg has been selected:
        function determinePeg(app)
            
            % Get click coordinates:
            mouseCoord = get(gca,'CurrentPoint');
            
            % Compute distance from each colored peg centre:
            deltaX = repmat(app.pegX - mouseCoord(1,1),1,2);
            deltaY = repmat(app.pegY - mouseCoord(1,2),1,3)';
            distanceMap = sqrt(deltaX.^2+deltaY.^2);
            minDistance = min(min(distanceMap));
            app.colorIdx = [];
            
            % Check if click happened within peg boundaries:
            if minDistance <= app.pegRadius
                % Determine closest colored peg:
                app.colorIdx = find(distanceMap==minDistance);
            end
            
            % Check if peg has been selected:
            if isempty(app.colorIdx)
                % Don't do nothing:
                set(gcf,'WindowButtonMotionFcn',@app.emptyFcn);
                set(gcf,'WindowButtonUpFcn',@app.emptyFcn);
            else
                % Store peg position and increase visibility:
                app.defaultPegPosition = [app.coloredPeg(app.colorIdx).obj.XData;app.coloredPeg(app.colorIdx).obj.YData];
                uistack(app.coloredPeg(app.colorIdx).obj,'top');
                
                % Setup functions for drag and click release:
                set(gcf,'WindowButtonMotionFcn',@app.movePeg);
                set(gcf,'WindowButtonUpFcn',@app.makeMove);
            end
        end
        
        %% Make the colored peg follow your cursor:
        function movePeg(app,~,~)
            
            % Get click coordinates:
            mouseCoord = get(gca,'CurrentPoint');
            
            % Update X and Y coordiantes of the selected colored peg:
            app.coloredPeg(app.colorIdx).obj.XData = [mouseCoord(1,1)-app.pegRadius mouseCoord(1,1)+app.pegRadius];
            app.coloredPeg(app.colorIdx).obj.YData = [mouseCoord(1,2)-app.pegRadius mouseCoord(1,2)+app.pegRadius];
        end
        
        %% Place a colored peg on the board:
        function makeMove(app,~,~)
            
            % Make the peg stop following cursor and place it back:
            set(gcf,'WindowButtonMotionFcn',@app.emptyFcn);
            app.coloredPeg(app.colorIdx).obj.XData = app.defaultPegPosition(1,:);
            app.coloredPeg(app.colorIdx).obj.YData = app.defaultPegPosition(2,:);
            
            % Get click coordinates:
            mouseCoord = get(gca,'CurrentPoint');
            
            % Compute distance from each peg hole centre:
            deltaX = app.pegHoleX - mouseCoord(1,1);
            deltaY = repmat(app.pegHoleY(app.currentRow) - mouseCoord(1,2),1,4)';
            distanceMap = sqrt(deltaX.^2+deltaY.^2);
            minDistance = min(min(distanceMap));
            
            % Check if click happened within peg boundaries:
            if minDistance <= app.pegRadius
                % Determine closest peg hole and update its color:
                app.holePeg(app.currentRow,find(distanceMap==minDistance)).obj.CData = app.pegGRX(app.colorIdx).CData;
                
                % Update current solution attempt:
                app.attempt(find(distanceMap==minDistance)) = app.colorIdx;
                
                % Check if all four pegs have been placed:
                checkEndTurn(app);
            end
        end
        
        %% Check if all four pegs have been placed:
        function checkEndTurn(app)
            
            % Check if the current attempt has been completed:
            if ~any(app.attempt == length(app.pegColors)+1)
                % Check if the attempt is winning:
                if isequal(app.attempt,app.solution)
                    % Stop game and show end-game winning splashStart art
                    app.isPlaying = 0;
                    showResult(app,'YouWin');
                % Checks if the limit number of tries has been exceeded:
                elseif app.currentRow == app.nRows
                    % Stop game and show end-game losing splashStart art
                    app.isPlaying = 0;
                    showResult(app,'GameOver');
                end
                
                % Returns the feedback pegs:
                getFeedback(app);
                
                % Update engine data:
                app.allSolutions(all(app.allSolutions==app.attempt,2),:) = [];
                
                % Update game parameters:
                app.currentRow = app.currentRow + 1;
                app.attempt = 7*ones(1,4);
            end
        end
        
        %% Compute the feedback pegs (white/black) based on attempt and solution:
        function getFeedback(app)
            
            % Count per color of the attempt code:
            attemptCount = countSolutionNumber(app,app.attempt);
            
            % Evaluation of the number of white and black pegs:
            whiteDotsAll = sum(min([attemptCount(app.solutionCount>0);app.solutionCount(app.solutionCount>0)]));
            blackDots = length(find(app.solution-app.attempt==0));
            whiteDots = whiteDotsAll-blackDots;
            
            % Update result on GUI:
            idx = 0;
            for i = 1 : blackDots
                % Color first pegs as black:
                idx = idx + 1;
                app.resultsPeg(app.currentRow,idx).obj.FaceColor = 'k';
            end
            for i = 1 : whiteDots
                % Color following pegs as white:
                idx = idx + 1;
                app.resultsPeg(app.currentRow,idx).obj.FaceColor = 'w';
            end
            
            % Update engine list of possible solutions:
            excludeIdx = excludeImpossibleSolutions(app,app.possibleSolutions,app.attempt,[blackDots,whiteDots]);
            app.possibleSolutions(excludeIdx,:) = [];
        end
        
        %% Get a suggestion from the engine for the next move:
        function useEngine(app)
            
            if app.currentRow == 1
                % Default hint at first try:
                optimalGuess = [1 1 2 2];
            elseif isempty(app.currentRow)
                % Game is not started yet:
                return
            elseif size(app.possibleSolutions,1) == 1
                % Returns the solution when the number of possible solutions is equal to one:
                optimalGuess = app.possibleSolutions;
            else
                % Use Donald Knuth algorithm for next attempt hint:
                optimalGuess = getOptimalGuess(app);
            end
            
            % Represent the suggestion on the board:
            for j = 1 :4
                app.holePeg(app.currentRow,j).obj.CData = .5*app.pegGRX(optimalGuess(j)).CData;
            end
        end
        
        %% Determine which button has been selected:
        function makeAction(app)
            
            mouseCoord = get(gca,'CurrentPoint');
            if mouseCoord(1,2) > app.button1Position(2) && mouseCoord(1,2)<app.button1Position(2)+app.button1Position(4)
                % Change temporarily Exit button's color if pressed and redirects:
                createButton(app,1,2);
                pause(.1);
                close(app)
            elseif mouseCoord(1,2) > app.button2Position(2) && mouseCoord(1,2)<app.button2Position(2)+app.button2Position(4) && app.isPlaying
                % Change temporarily Hint button's color if pressed and redirects:
                createButton(app,2,2);
                pause(.1);
                useEngine(app)
                createButton(app,2,1);
            elseif mouseCoord(1,2) > app.button3Position(2) && mouseCoord(1,2)<app.button3Position(2)+app.button3Position(4)
                % Change temporarily Play button's color if pressed and redirects:
                createButton(app,3,2);
                pause(.1);
                newGame(app)
                createButton(app,3,1);
            end
        end
        
        %% Exclude impossible solutions based on the attempt feedbacks:
        function excludeIdx = excludeImpossibleSolutions(app,possibleSolutions,guess,feedback)
            
            % Initialise the indeces corresponding to impossible solutions:
            excludeIdx = [];
            
            for j = 1 : size(possibleSolutions,1)
                
                % Assesses the feedback for every possible solution, based on the current attempt on the board:
                % [See getFeedback(app) for more informations]
                guessCount = countSolutionNumber(app,guess);
                currentSolutionCount = countSolutionNumber(app,possibleSolutions(j,:));
                whiteDotsAll = sum(min([guessCount(currentSolutionCount>0);currentSolutionCount(currentSolutionCount>0)]));
                blackDots = length(find(possibleSolutions(j,:)-guess==0));
                whiteDots = whiteDotsAll-blackDots;
                
                % Excludes the codes that have a different feedback from
                % the one on the board:
                if ~isequal([blackDots,whiteDots],feedback)
                    excludeIdx(end+1) = j;
                end
            end
        end
        
        %% Computes the engine hint code as the one reducing most the possible solutions:
        function optimalGuess = getOptimalGuess(app)
            
            % Initialise variables:
            N = size(app.allSolutions,1);
            K = size(app.feedbackList,1);
            guessList = zeros(N,app.nCols+1);
            
            % Start loop over every possible non-attempted code (minmax):
            for j = 1 : N
                result = zeros(K,1);
                guess = app.allSolutions(j,:);
                for k = 1 : K
                    % Computes the number of exclusion for every possible feedback:
                    feedback = app.feedbackList(k,:);
                    excludeIdx = excludeImpossibleSolutions(app,app.possibleSolutions,guess,feedback);
                    result(k) = length(excludeIdx);
                end
                % Retains the lowest value of exclusion:
                guessList(j,:) = [guess,min(result)];
            end
            
            % Use the code with more exclusions as next hint:
            guessList = sortrows(guessList,app.nCols+1);
            optimalGuess = guessList(end,1:4);
        end
        
        %% Computes all the possible 1296 solutions:
        function initialCombinations = getAllCombinations(app)
            
            % Permutation with repetitions:
            x = 1:length(app.pegColors);
            C = cell(app.nCols,1);
            [C{:}] = ndgrid(x);
            y = cellfun(@(x){x(:)},C);
            initialCombinations = [y{:}];
        end
        
        %% Shows end-game splashStart art for win and loss:
        function showResult(app,flag)
            
            % Position parameters for splashStart art:
            xImg = [app.msgPosition(1)+app.resultBaseGap,app.msgPosition(1)+app.msgPosition(3)-app.resultBaseGap];
            yImg = [app.msgPosition(2)+app.msgPosition(4)-app.resultBaseGap,app.msgPosition(2)+app.resultLowerGap];
            
            % Position parameters for solution:
            onePegSpace = 2*app.pegRadius+app.pegColGap;
            xSolutionGap = (app.msgPosition(3)-(app.nCols-1)*(onePegSpace))/2;
            xSolution = app.msgPosition(1)+xSolutionGap+(0:3)*onePegSpace;
            ySolution = app.msgPosition(2)+app.resultLowerGap/2;
            
            % splashStart art:
            app.splashEnd.box = rectangle(app.gameAxis,'Position',app.msgPosition,'Curvature',0.3,'FaceColor','w','EdgeColor','k','LineWidth',10);
            if app.isData
                app.splashEnd.text = image(app.gameAxis,'XData',xImg,'YData',yImg,'CData',app.GRXdata.(flag));
                colormap(gray)
            end
            
            % Colored peg solution code on GUI:
            for j = 1 : app.nCols
                app.splashEnd.peg(j).obj = image(app.gameAxis,'XData',xSolution(j)+app.pegRadius*[-1 1],'YData',ySolution+app.pegRadius*[-1 1],'CData',app.pegGRX(app.solution(j)).CData,'AlphaData',app.alphamap,'CDataMapping','scaled','AlphaDataMapping','scaled','Interpolation','bilinear');
            end
            
            % Detailed borders:
            th = 0:.01:pi;
            L = length(th);
            Xrow = [xSolution(1)+(app.pegRadius+app.pegRowGap/2)*cos(th+pi/2), xSolution(end)+(app.pegRadius+app.pegRowGap/2)*cos(th-pi/2)]; Xrow = [Xrow,Xrow(1)];
            Yrow = [ySolution+(app.pegRadius+app.pegRowGap/2)*sin(th+pi/2), ySolution+(app.pegRadius+app.pegRowGap/2)*sin(th-pi/2)]; Yrow = [Yrow,Yrow(1)];
            app.splashEnd.details(1).line = plot(app.gameAxis,Xrow(L:L+1),Yrow(L:L+1)+2,'Color',.5*[0.8314 0.6863 0.2157],'LineWidth',app.detailLineWidth);
            app.splashEnd.details(2).line = plot(app.gameAxis,Xrow(L:L+1),Yrow(L:L+1)-2,'Color',.5*[0.8314 0.6863 0.2157],'LineWidth',app.detailLineWidth);
            app.splashEnd.details(3).line = plot(app.gameAxis,Xrow,Yrow,'-o','Color',[0.8314 0.6863 0.2157],'LineWidth',2,'MarkerSize',app.detailLineWidth+1);
            app.splashEnd.details(4).line = plot(app.gameAxis,Xrow(L:L+1),Yrow(2*L:2*L+1)+2,'Color',.5*[0.8314 0.6863 0.2157],'LineWidth',app.detailLineWidth);
            app.splashEnd.details(5).line = plot(app.gameAxis,Xrow(L:L+1),Yrow(2*L:2*L+1)-2,'Color',.5*[0.8314 0.6863 0.2157],'LineWidth',app.detailLineWidth);
            
            % Update game parameter:
            app.multipleGame = 1;
        end
        
        % =================================================================
        % GUI FUNCTIONS
        % =================================================================
        
        %% Create the GUI and initialise components and parameters:
        function initialiseGUI(app)
            
            % Panels position parameters [px]:
            leftWidth = 2*app.borderWidth + (app.nCols+1)*app.pegColGap + 2*app.nCols*app.pegRadius;
            leftBottomHeight = 2*app.borderWidth + 2*app.initialGap + 3*app.pegRowGap + 2*2*app.pegRadius;
            leftUpHeight = 2*app.borderWidth + 2*app.initialGap + (app.nRows+1)*app.pegRowGap + 2*app.nRows*app.pegRadius;
            rightWidth = 2*app.borderWidth + 2*app.resultsShift + 2*app.pegColGap;
            buttonHeight = leftBottomHeight/3;
            msgHorGap = app.borderWidth + app.pegColGap;
            msgVerGap = (leftBottomHeight + leftUpHeight - (leftWidth + rightWidth - 2*msgHorGap))/2;
            
            % Panel position vectors [px]:
            leftBottomPosition = [0,0,leftWidth,leftBottomHeight];
            leftUpPosition = [0,leftBottomHeight,leftWidth,leftUpHeight];
            rightPosition = [leftWidth,leftBottomHeight,rightWidth,leftUpHeight];
            figurePosition = [500,50,leftWidth+rightWidth,leftBottomHeight+leftUpHeight];
            app.button1Position = [leftWidth,0,rightWidth,buttonHeight];
            app.button2Position = [leftWidth,buttonHeight,rightWidth,buttonHeight];
            app.button3Position = [leftWidth,2*buttonHeight,rightWidth,buttonHeight];
            app.msgPosition = [msgHorGap,msgVerGap,figurePosition(3)-2*msgHorGap,figurePosition(4)-2*msgVerGap];
            
            % Main figure:
            mainWindow = figure('Visible','off');
            mainWindow.Position = figurePosition;
            mainWindow.Name = 'Mastermind';
            mainWindow.MenuBar = 'None';
            mainWindow.NumberTitle = 'off';
            mainWindow.Renderer = 'zbuffer';
            mainWindow.Resize = 'off';
            mainWindow.Color = app.windowColor;
            mainWindow.Visible = 'on';
            set(gcf,'CloseRequestFcn',@app.close)
            
            % Axis:
            app.gameAxis = axes(mainWindow);
            app.gameAxis.Position = [0 0 1 1];
            hold(app.gameAxis,'all')
            axis(app.gameAxis,'off');
            xlim(app.gameAxis,[0 leftWidth+rightWidth]);
            ylim(app.gameAxis,[0 leftBottomHeight+leftUpHeight]);
            disableDefaultInteractivity(app.gameAxis);
            
            % Drawing panels borders:
            % [Credits: Brian Moore - https://it.mathworks.com/matlabcentral/profile/authors/3690384]
            [Xlb,Ylb] = getPatchCoord(app,leftBottomPosition);
            [Xlu,Ylu] = getPatchCoord(app,leftUpPosition);
            [Xru,Yru] = getPatchCoord(app,rightPosition);
            patchColor = permute([app.darkColor',app.lightColor'],[3 2 1]);
            patch(app.gameAxis,Xlb,Ylb,patchColor,'EdgeColor','None');
            patch(app.gameAxis,Xlu,Ylu,patchColor,'EdgeColor','None');
            patch(app.gameAxis,Xru,Yru,patchColor,'EdgeColor','None');
            
            % Create buttons:
            createButton(app,1,1);
            createButton(app,2,1);
            createButton(app,3,1);
            
            % Add graphic elements:
            pegGraphicsData(app);
            getHoleGrid(app,leftUpPosition);
            getPegGrid(app,leftBottomPosition);
            getResultsGrid(app,rightPosition);
            plotResultsPegs(app);
            plotMainPegs(app);
            
            % Check if data.mat is found:
            gamePath = fileparts(matlab.desktop.editor.getActiveFilename);
            app.isData = any(contains({dir(gamePath).name},'data.mat'));
            
            % Start-game splashStart art:
            if app.isData
                app.GRXdata = load(fullfile(gamePath,'data.mat'));
                leftX = app.borderWidth;
                rightX = leftWidth - leftX;
                lowY = leftBottomHeight + app.borderWidth + app.initialGap;
                upY = figurePosition(4) - app.borderWidth - app.initialGap;
                app.splashStart = image(app.gameAxis,'XData',[leftX,rightX],'YData',[upY,lowY],'CData',app.GRXdata.CData,'AlphaData',app.GRXdata.alphamap);
            end
            
            % Activate GUI:
            mainWindow.Visible = 'on';
            set(gcf,'WindowButtonDownFcn',@app.determineAction)
            app.isPlaying = 0;
        end
        
        %% Compute coordinates of the peg holes:
        function getHoleGrid(app,position)
            
            % Distance between holes:
            pegDistanceX = 2*app.pegRadius + app.pegColGap;
            pegDistanceY = 2*app.pegRadius + app.pegRowGap;
            
            % Peg hole initial position:
            X = (position(1) + app.borderWidth + app.pegColGap + app.pegRadius)*ones(app.nCols,1);
            Y = (position(2) + app.borderWidth + app.initialGap + app.pegRowGap + app.pegRadius)*ones(app.nRows,1);
            
            % Iterate over the grid layout:
            for j = 2 : app.nRows
                if j <= app.nCols
                    X(j:end) = X(j:end) + pegDistanceX;
                end
                Y(j:end) = Y(j:end) + pegDistanceY;
            end
            
            % Set coordinates as properties:
            app.pegHoleX = X;
            app.pegHoleY = Y;
        end
        
        %% Compute coordinates of the colored pegs:
        function getPegGrid(app,position)
            
            % Distance between pegs:
            numRow = 2;
            numCol = length(app.pegColors)/numRow;
            pegDistanceX = (position(3)-2*app.borderWidth+2*app.pegRadius)/(numCol+1);
            pegDistanceY = (position(4)-2*app.borderWidth+3*app.pegRadius/2)/(numRow+1);
            
            % Peg hole initial position:
            X = (app.borderWidth + pegDistanceX-app.pegRadius)*ones(numCol,1);
            Y = (app.borderWidth + pegDistanceY-3*app.pegRadius/4)*ones(numRow,1);
            
            % Iterate over the grid layout:
            for j = 2 : numCol
                if j <= numRow
                Y(j:end) = Y(j:end) + pegDistanceY;
                end
                X(j:end) = X(j:end) + pegDistanceX;
            end
            
            % Set coordinates as properties:
            app.pegX = X;
            app.pegY = Y;
        end
        
        %% Compute X-coordinate of the feedback pegs:
        function getResultsGrid(app,position)
            X = position(1) + position(3)/2;
            app.resultsX = X;
        end
        
        %% Computes CData and AlphaMap of shadowed pegs:
        function pegGraphicsData(app)
            
            % Initialise parameters:
            r = 10*app.pegRadius;
            allColors = [app.pegColors,'#ffffff'];
            sign = -1;
            
            % Number of pixels along one dimension:
            n = 2*r+1;
            
            % Iterate over all colors plus gray for peg holes:
            for i = 1 : length(allColors)
                
                % Last index refers to peg holes:
                if i == length(allColors)
                    sign = -sign;
                end
                
                % Converts hex string to RGB normalised:
                color = allColors{i}(2:end);
                RGB = [hex2dec(color(1:2))/255,hex2dec(color(3:4))/255,hex2dec(color(5:6))/255];
                
                % Initialise variables:
                CData = ones(n);
                AData = ones(n);
                
                % Assess the point of maximum luminosity:
                nHalf = ceil(n/2);
                lateralDisp = sign*ceil(nHalf*app.lightCoeffH);
                verticalDisp = sign*ceil(nHalf*app.lightCoeffV);
                centralSquare = [nHalf - lateralDisp , nHalf - verticalDisp];
                
                % Create CData and AlphaMap based on the distance from the point of maximum luminosity:
                for j = 1 : n
                    for k = 1 : n
                        % Check if the point lies within the peg:
                        distFromCenter = sqrt((j-nHalf)^2 + (k-nHalf)^2);
                        if distFromCenter <= r
                            % Assess CData values:
                            distFromLight = .92*sqrt((j-centralSquare(1))^2 + (k-centralSquare(2))^2);
                            idx = min(1,1 - distFromLight/r);
                            CData(j,k) = idx;
                        else
                            % Set transparency for points outside the peg:
                            AData(j,k) = 0;
                        end
                    end
                end
                
                % Delete corners:
                AData(nHalf,nHalf+r) = 0;
                AData(nHalf,nHalf-r) = 0;
                AData(nHalf+r,nHalf) = 0;
                AData(nHalf-r,nHalf) = 0;
                
                % Add colors to CData:
                CData = repmat(CData,1,1,3);
                CData(:,:,1) = CData(:,:,2)*RGB(1);
                CData(:,:,2) = CData(:,:,2)*RGB(2);
                CData(:,:,3) = CData(:,:,3)*RGB(3);
                
                % Store as properties:
                app.pegGRX(i).CData = CData;
                app.alphamap = AData;
            end
        end
        
        %% Computes the coordinates of the borders for the patch:
        function [X,Y] = getPatchCoord(app,position)
            
            % Corner positions:
            x0 = position(1);
            y0 = position(2);
            xN = position(1) + position(3);
            yN = position(2) + position(4);
            
            % Border Width:
            bw = app.borderWidth;
            
            % Output coordinates:
            X = [x0,x0,xN,xN-bw,x0+bw,x0+bw;x0,xN,xN,xN-bw,xN-bw,x0+bw]';
            Y = [y0,yN,yN,yN-bw,yN-bw,y0+bw;y0,y0,yN,yN-bw,y0+bw,y0+bw]';
        end
        
        %% Represent peg holes on the GUI and add some details:
        function plotPegHoles(app)
            
            % Initialise parameters:
            th = 0:.01:pi;
            L = length(th);
            
            % Iterate the plot over the holes grid:
            for i = 1 : app.nRows
                y = app.pegHoleY(i) + app.pegRadius*[-1 1];
                for j = 1 : app.nCols
                    x = app.pegHoleX(j) + app.pegRadius*[-1 1];
                    app.holePeg(i,j).obj = image(app.gameAxis,'XData',x,'YData',y,'CData',app.pegGRX(7).CData,'AlphaData',app.alphamap,'CDataMapping','scaled','AlphaDataMapping','scaled','Interpolation','bilinear');
                end
                
                % Add details around peg holes:
                Xrow = [app.pegHoleX(1)+(app.pegRadius+app.pegRowGap/2)*cos(th+pi/2), app.pegHoleX(end)+(app.pegRadius+app.pegRowGap/2)*cos(th-pi/2)]; Xrow = [Xrow,Xrow(1)];
                Yrow = [app.pegHoleY(i)+(app.pegRadius+app.pegRowGap/2)*sin(th+pi/2), app.pegHoleY(i)+(app.pegRadius+app.pegRowGap/2)*sin(th-pi/2)]; Yrow = [Yrow,Yrow(1)];
                plot(app.gameAxis,Xrow(L:L+1),Yrow(L:L+1)+4,'Color',app.lightColor,'LineWidth',app.detailLineWidth)
                plot(app.gameAxis,Xrow(L:L+1),Yrow(L:L+1)-4,'Color',app.lightColor,'LineWidth',app.detailLineWidth)
                plot(app.gameAxis,Xrow,Yrow,'-o','Color',app.darkColor,'LineWidth',2,'MarkerSize',app.detailLineWidth+1)
            end
            
            % Conclude details with upper horizontal line:
            plot(app.gameAxis,Xrow(L:L+1),Yrow(2*L:2*L+1)+4,'Color',app.lightColor,'LineWidth',app.detailLineWidth)
            plot(app.gameAxis,Xrow(L:L+1),Yrow(2*L:2*L+1)-4,'Color',app.lightColor,'LineWidth',app.detailLineWidth)
        end
        
        %% Represent colored pegs on the GUI:
        function plotMainPegs(app)
            
            % Initialise parameters for pegs grid position:
            numRow = 2;
            numCol = length(app.pegColors)/numRow;
            pegID = 0;
            
            % Iterate the plot over the grid:
            for i = 1 : numRow
                    y = app.pegY(i) + app.pegRadius*[-1 1];
                for j = 1 : numCol
                    pegID = pegID+1;
                    x = app.pegX(j) + app.pegRadius*[-1 1];
                    app.coloredPeg(pegID).obj = image(app.gameAxis,'XData',x,'YData',y,'CData',app.pegGRX(pegID).CData,'AlphaData',app.alphamap,'CDataMapping','scaled','AlphaDataMapping','scaled','Interpolation','bilinear');
                end
            end
        end
        
        %% Represent the feedback pegs on the GUI and add some details:
        function plotResultsPegs(app)
            
            % Position shift of each peg from middle line of the panel:
            positionShift = (app.resultsShift)*ones(1,2).*[(dec2bin(0:3)-'0')*2-1];
            
            % Iterate the plot over the grid:
            for i = 1 : app.nRows
                for j = 1 : app.nCols
                    x = app.resultsX + positionShift(j,1)-app.resultsRadius;
                    y = app.pegHoleY(i) + positionShift(j,2)-app.resultsRadius;
                    app.resultsPeg(i,j).obj = rectangle('Position',[x,y,2*app.resultsRadius,2*app.resultsRadius],'Curvature',[1 1],'FaceColor',[.3 .3 .3],'EdgeColor',[.3 .3 .3]);
                end
                
                % Add details:
                rectangle('Position',[app.resultsX-2*app.resultsShift-2,app.pegHoleY(i)-2*app.resultsShift-2,4*app.resultsShift+4,4*app.resultsShift+4],'Curvature',.3,'EdgeColor','w','linewidth',3.5);
                rectangle('Position',[app.resultsX-2*app.resultsShift-2,app.pegHoleY(i)-2*app.resultsShift-2,4*app.resultsShift+4,4*app.resultsShift+4],'Curvature',.3,'EdgeColor',app.windowColor,'linewidth',1.75);
            end
        end
        
        %% Represent the buttons on the GUI:
        function createButton(app,buttonIdx,coeff)
            
            % Initialise position, color and label of each button:
            switch buttonIdx
                case 1
                    position = app.button1Position;
                    c(:,:,1) = 1/coeff; c(:,:,2) = 0; c(:,:,3) = 0;
                    str = 'EXIT'; horSpaceCoeff = 1.5;
                case 2
                    position = app.button2Position;
                    c(:,:,1) = 0; c(:,:,2) = 1/coeff; c(:,:,3) = 1/coeff;
                    str = 'HINT'; horSpaceCoeff = 1.1;
                case 3
                    position = app.button3Position;
                    c(:,:,1) = 0; c(:,:,2) = 1/coeff; c(:,:,3) = 0;
                    str = 'PLAY'; horSpaceCoeff = 1.4;
            end
            
            % Button position:
            Xb = [position(1),position(1)+position(3),position(1)+position(3),position(1)];
            Yb = [position(2),position(2),position(2)+position(4),position(2)+position(4)];
            
            % Get coordinates and colors for the button borders:
            % [Credits: Brian Moore - https://it.mathworks.com/matlabcentral/profile/authors/3690384]
            [Xb1,Yb1] = getPatchCoord(app,position);
            patchColor = permute([app.lightColor',app.darkColor'],[3 2 1]);
            
            % GUI elements:
            patch(app.gameAxis,Xb,Yb,app.windowColor.*permute(c,[2,3,1]),'EdgeColor','None');
            patch(app.gameAxis,Xb1,Yb1,patchColor.*c,'EdgeColor','None');
            text(position(1)+horSpaceCoeff*app.borderWidth,position(2)+2.5*app.borderWidth,str,'Color',abs(1-coeff)*[1 1 1],'FontSize',18,'FontName','Arial','FontWeight','bold','FontName','Ink Free')
        end
        
        % =================================================================
        % MISCELLANEOUS
        % =================================================================
        
        %% Close the GUI:
        function close(app,~,~)
            delete(app);
            delete(gcf);
        end
        
        %% Do nothing:
        function emptyFcn(app,~,~)
        end
        
    end
    
   methods (Access = public)
       
       %% Start application:
       function app = Mastermind
           initialiseGUI(app)
       end
       
   end
end