#include "firesim/FireSim.h"
#include "fluidsim/FluidSim.h"
#include "baseinterface/SimulationSettings.h"

#include <chrono>
#include <SFML/Graphics.hpp>
#include <iostream>

void makeImageFromArray(std::vector<std::vector<SIM::colour >>& pixelArray,int size, sf::Image &simImage, int windowSize) {
    int scale = windowSize / size;
    for (unsigned int i = 0; i < size - 1; ++i) {
        for (unsigned int j = 0; j < size - 1; ++j) {
            auto& pix = pixelArray[j][i];
            auto pixcol = sf::Color{ pix.r,pix.g,pix.b };
            for (int is = 0; is < scale; ++is)
            {
            	for (int js = 0; js < scale; ++js)
	            {
                    simImage.setPixel(i*scale+is, j*scale+js, pixcol);
	            }
            }
            

        }
    }
}


int main() {
    int size = 512;
    int scale = 2;
    int windowSize = size*scale;
    sf::Image simImage;
    simImage.create(windowSize, windowSize);
    sf::Clock gameclock;
    sf::RenderWindow window(sf::VideoMode(windowSize, windowSize), "simulations");
    auto settings = std::make_unique<SIM::SimulationSettings>(size, false);
    settings->setDiffusion(0.000001);
    settings->setViscosity(0.000001);
    SIM::FluidSim sim{ std::move(settings) };

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            switch (event.type)
            {
            case sf::Event::Closed:
                window.close();
                break;
            case  sf::Event::MouseButtonPressed:
                if (event.mouseButton.button == sf::Mouse::Right)
                {
                    sim.handleClick(false, event.mouseButton.x/scale, event.mouseButton.y/scale);
                }else
                {
                    sim.handleClick(true, event.mouseButton.x / scale, event.mouseButton.y / scale);
                }
                break;
            }
        }
        //if (gameclock.getElapsedTime().asSeconds() > 0.013) {
            //sim.advance(gameclock.getElapsedTime().asSeconds());
            
            sim.advance(0.013);
            window.clear();

            makeImageFromArray(sim.getCurrentState(), size, simImage, windowSize);
            sf::Texture texture;
            texture.loadFromImage(simImage);
            sf::Sprite sprite;
            sprite.setTexture(texture, true);
            window.draw(sprite);
            std::cout << gameclock.getElapsedTime().asSeconds() << std::endl;
            gameclock.restart();
            window.display();
        //}

    }

    return 0;
    
}
