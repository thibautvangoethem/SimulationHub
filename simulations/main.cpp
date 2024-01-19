#include "core/baseinterface/SimulationSettings.h"
#include "sfmlFrontend/SFMLSimInterface.h"
#include "sfmlFrontend/FireSimImpl.h"
#include "sfmlFrontend/FluidSimImpl.h"
#include "sfmlFrontend/RayTracerImpl.h"
#include "sfmlFrontend/SlimeSimimpl.h"

#include <chrono>
#include <functional>
#include <SFML/Graphics.hpp>
#include <iostream>

void runSimulation(sf::RenderWindow& window, int scale, int size, int windowSize, std::unique_ptr<FSIM::SFMLSimInterface> sim)
{
	auto currentTimeCount = 0.0;
	auto frameCounter = 0;

	bool timeBased = sim->getInternalSettings()->isTimeBased();
	bool unlocked = sim->getInternalSettings()->isUnlockedTime();
	double tickrate = sim->getInternalSettings()->getConstantTickrate();
	sf::Image simImage;
	simImage.create(windowSize, windowSize);
	sf::Clock gameclock;
	while (window.isOpen())
	{
		if (unlocked || gameclock.getElapsedTime().asSeconds() > 0.13) {
			sf::Event event;
			while (window.pollEvent(event))
			{
				switch (event.type)
				{
				case sf::Event::Closed:
					window.close();
					break;
				case sf::Event::MouseButtonPressed:
					if (event.mouseButton.button == sf::Mouse::Right)
					{
						sim->click(false, event.mouseButton.x / scale, event.mouseButton.y / scale);
					}
					else
					{
						sim->click(true, event.mouseButton.x / scale, event.mouseButton.y / scale);
					}
					break;
				case sf::Event::KeyPressed:
					if (event.key.code == sf::Keyboard::Escape)
						return;
					break;
				}
			}
			sim->advanceSim(timeBased ? gameclock.getElapsedTime().asSeconds() : tickrate);
			gameclock.restart();
			sim->updateImage(simImage);
			window.clear();
			//makeImageFromArray(sim->getCurrentState(), size, simImage, windowSize);
			sf::Texture texture;
			texture.loadFromImage(simImage);
			sf::Sprite sprite;
			sprite.setTexture(texture, true);
			window.draw(sprite);
			frameCounter += 1;
			currentTimeCount += gameclock.getElapsedTime().asSeconds();
			if (frameCounter == 1)
			{
				std::cout << "first frame: " << currentTimeCount << "s" << std::endl;
			}
			if(frameCounter==500)
			{
				std::cout <<"500 frame average: " << currentTimeCount/500 <<"s"<< std::endl;
				currentTimeCount = 0;
				frameCounter = 1;
			}

			
			window.display();
		}
	}
}

int main()
{
	//Dont feel like dealing with adding resources to an executable atm
	std::string fontLocation = R"(C:\Users\thiba\source\repos\FluidSim\arial\arial.ttf)";
	int size = 256;
	int scale = 4;
	int windowSize = size * scale;
	sf::RenderWindow window(sf::VideoMode(windowSize, windowSize), "simulations");
	auto settings = std::make_shared<SIM::SimulationSettings>(size, false);
	settings->setConstantTickrate(0.01);
	settings->setunlockedTime(true);
	settings->setDiffusion(0.000001);
	settings->setViscosity(0.000001);
	settings->setTimeBased(false);


	std::vector<std::function<std::unique_ptr<FSIM::SFMLSimInterface>()>> choices;

	//some might say there are better ways to do this, I agree with them. Luckily this will not be a long term project right?(written 2023/04/02).
	std::function<std::unique_ptr<FSIM::SFMLSimInterface>()> func1 = std::function(
		[&]() { return std::make_unique<FSIM::FireSimImpl>(FSIM::FireSimImpl{ settings,scale }); });
	choices.emplace_back(func1);
	std::function<std::unique_ptr<FSIM::SFMLSimInterface>()> func2 = std::function(
		[&]() { return std::make_unique<FSIM::FluidSimImpl>(FSIM::FluidSimImpl{ settings,scale }); });
	choices.emplace_back(func2);
	std::function<std::unique_ptr<FSIM::SFMLSimInterface>()> func3 = std::function(
		[&]() { return std::make_unique<FSIM::RayTracerImpl>(FSIM::RayTracerImpl{ settings ,scale}); });
	choices.emplace_back(func3);
	std::function<std::unique_ptr<FSIM::SFMLSimInterface>()> func4 = std::function(
		[&]() { return std::make_unique<FSIM::SlimeSimImpl>(FSIM::SlimeSimImpl{ settings ,scale }); });
	choices.emplace_back(func4);

	sf::Font font;
	font.loadFromFile(fontLocation);
	int split = windowSize / 4;
	int implemented = choices.size();
	bool running = false;
	int simIndexRunning=0;
	while (window.isOpen())
	{
		window.clear();
		sf::Event event;
		while (window.pollEvent(event))
		{
			switch (event.type)
			{
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::MouseButtonPressed:
				
				simIndexRunning = event.mouseButton.x / split + (event.mouseButton.y / split) * 4;
				if (simIndexRunning < implemented) {
					running = true;
					std::cout << "Running simulation: " << simIndexRunning << std::endl;
				}
				else
				{
					std::cout << "This index is not implemented yet, maybe 16 is a bit too ambitious here.";
				}
				break;
			}
		}
		if (running) {
			runSimulation(window, scale, size, windowSize, choices[simIndexRunning]());
			running = false;
		}
		else {
			for (auto i = 0; i < 4; ++i)
			{
				for (auto j = 0; j < 4; ++j)
				{
					int index = i + j * 4;
					sf::RectangleShape rectangle;
					rectangle.setSize(sf::Vector2f(split, split));
					rectangle.setOutlineColor(sf::Color::Blue);
					rectangle.setOutlineThickness(5);
					rectangle.setPosition(i * split, j * split);
					window.draw(rectangle);
					if (index < implemented) {
						sf::Text text(std::to_string(index), font);
						text.setCharacterSize(30);
						text.setStyle(sf::Text::Bold);
						text.setFillColor(sf::Color::Black);
						text.setPosition(i * split, j * split);
						window.draw(text);
					}
					else
					{
						sf::Vertex line1[2];
						line1[0].position = sf::Vector2f(i * split, j * split);
						line1[0].color = sf::Color::Red;
						line1[1].position = sf::Vector2f((i + 1) * split, (j + 1) * split);
						line1[1].color = sf::Color::Red;

						sf::Vertex line2[2];
						line2[0].position = sf::Vector2f((i + 1) * split, j * split);
						line2[0].color = sf::Color::Red;
						line2[1].position = sf::Vector2f(i * split, (j + 1) * split);
						line2[1].color = sf::Color::Red;

						window.draw(line1, 2, sf::Lines);
						window.draw(line2, 2, sf::Lines);
					}
				}
			}
			window.display();
		}
	}
	return 0;
}

