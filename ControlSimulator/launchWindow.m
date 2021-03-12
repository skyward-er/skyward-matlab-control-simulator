classdef launchWindow < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        launchButton                 matlab.ui.control.Button
        PressthebuttontolaunchLabel  matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: launchButton
        function launchButtonPushed(app, ~)
%             launch = true;
%             assignin('base', 'launchFlag', launch);
            fopen(strcat(pwd, '/launchFlag.txt'),'w+');
            
            delete(app);
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 382 149];
            app.UIFigure.Name = 'MATLAB App';

            % Create launchButton
            app.launchButton = uibutton(app.UIFigure, 'push');
            app.launchButton.ButtonPushedFcn = createCallbackFcn(app, @launchButtonPushed, true);
            app.launchButton.Position = [142 33 100 22];
            app.launchButton.Text = 'launch';

            % Create PressthebuttontolaunchLabel
            app.PressthebuttontolaunchLabel = uilabel(app.UIFigure);
            app.PressthebuttontolaunchLabel.Position = [119 93 145 22];
            app.PressthebuttontolaunchLabel.Text = 'Press the button to launch';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = launchWindow

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end