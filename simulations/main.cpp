#include "baseinterface/Simulation.h"
#include "firesim/FireSim.h"
#include "fluidsim/FluidSim.h"
#include "baseinterface/SimulationSettings.h"

#include <chrono>
#include <functional>
#include <SFML/Graphics.hpp>
#include <iostream>

void makeImageFromArray(std::vector<std::vector<SIM::colour>>& pixelArray, int size, sf::Image& simImage,
	int windowSize)
{
	const int scale = windowSize / size;
	for (unsigned int i = 0; i < size - 1; ++i)
	{
		for (unsigned int j = 0; j < size - 1; ++j)
		{
			const auto& pix = pixelArray[j][i];
			const auto pixcol = sf::Color{ pix.r, pix.g, pix.b };
			for (int is = 0; is < scale; ++is)
			{
				for (int js = 0; js < scale; ++js)
				{
					simImage.setPixel(i * scale + is, j * scale + js, pixcol);
				}
			}
		}
	}
}

void runSimulation(sf::RenderWindow& window, int scale, int size, int windowSize, std::unique_ptr<SIM::Simulation> sim)
{
	bool timeBased = sim->getSettings()->isTimeBased();
	bool unlocked = sim->getSettings()->isUnlockedTime();
	double tickrate = sim->getSettings()->getConstantTickrate();
	sf::Image simImage;
	simImage.create(windowSize, windowSize);
	sf::Clock gameclock;
	while (window.isOpen())
	{
		if (unlocked || gameclock.getElapsedTime().asSeconds() > tickrate) {
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
						sim->handleClick(false, event.mouseButton.x / scale, event.mouseButton.y / scale);
					}
					else
					{
						sim->handleClick(true, event.mouseButton.x / scale, event.mouseButton.y / scale);
					}
					break;
				case sf::Event::KeyPressed:
					if (event.key.code == sf::Keyboard::Escape)
						return;
					break;
				}
			}
			sim->advance(timeBased ? gameclock.getElapsedTime().asSeconds() : tickrate);
			window.clear();

			makeImageFromArray(sim->getCurrentState(), size, simImage, windowSize);
			sf::Texture texture;
			texture.loadFromImage(simImage);
			sf::Sprite sprite;
			sprite.setTexture(texture, true);
			window.draw(sprite);
			std::cout << gameclock.getElapsedTime().asSeconds() << std::endl;
			//Note: placing the restart here means that with a tickrate of 0.5s it will actually be a 0.5+simulation time refreshrate. however the tickrate of the simulation will still remain the chosen tickrate
			gameclock.restart();
			window.display();
		}
	}
}


int main()
{
	//Dont feel like dealing with adding resources to an executable atm
	std::string fontLocation = R"(C:\Users\thiba\source\repos\FluidSim\arial\arial.ttf)";
	int size = 512;
	int scale = 2;
	int windowSize = size * scale;
	sf::RenderWindow window(sf::VideoMode(windowSize, windowSize), "simulations");
	auto settings = std::make_shared<SIM::SimulationSettings>(size, false);
	settings->setConstantTickrate(0.05);
	settings->setunlockedTime(false);
	settings->setDiffusion(0.000001);
	settings->setViscosity(0.000001);


	std::vector<std::function<std::unique_ptr<SIM::Simulation>()>> choices;

	//some might say there are better ways to do this, I agree with them. Luckily this will not be a long term project right?(written 2023/04/02).
	std::function<std::unique_ptr<SIM::Simulation>()> func1 = std::function(
		[&]() { return std::make_unique<FSIM::FireSim>(FSIM::FireSim{ settings }); });
	choices.emplace_back(func1);
	std::function<std::unique_ptr<SIM::Simulation>()> func2 = std::function(
		[&]() { return std::make_unique<SIM::FluidSim>(SIM::FluidSim{ settings }); });
	choices.emplace_back(func2);
	sf::Font font;
	font.loadFromFile(fontLocation);
	int split = windowSize / 4;
	int implemented = choices.size();
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
				int index = event.mouseButton.x / split + (event.mouseButton.y / split) * 4;
				if (index < implemented) {
					std::cout << "Running simulation: " << index << std::endl;
					runSimulation(window, scale, size, windowSize, choices[index]());
				}
				else
				{
					std::cout << "This index is not implemented yet, maybe 16 is a bit too ambitious here.";
				}
				break;
			}
		}

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
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


	return 0;
}

