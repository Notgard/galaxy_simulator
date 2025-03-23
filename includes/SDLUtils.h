#pragma once

#include <SDL.h>
#include <SDL_ttf.h>

#include "Box.h"

class SDLUtils
{
public:
    SDLUtils() {};
    ~SDLUtils()
    {
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
    };

protected:
    bool graphical = false;
    bool render_tree_flag = false;
    bool render_center_of_mass_flag = false;
    bool render_pause_flag = false;

    Uint32 total_frame_ticks = 0;
    Uint32 ticks = 0;
    Uint64 perf_ticks = 0;

    SDL_Window *window = nullptr;
    SDL_Renderer *renderer = nullptr;
    TTF_Font *font;

    void draw_circle(int centerX, int centerY, SDL_Color color, float radius)
    {

        SDL_SetRenderDrawColor(this->renderer, color.r, color.g, color.b, color.a);

        // Loop through angles 0 to 2*PI, drawing points along the circle's circumference
        for (double angle = 0; angle < 2 * M_PI; angle += 0.05) // step size is small to smooth out the circle
        {
            int x = centerX + static_cast<int>(radius * cos(angle));
            int y = centerY + static_cast<int>(radius * sin(angle));
            SDL_RenderDrawPoint(this->renderer, x, y);
        }
    }

    void draw_box(const Box<float> &box, SDL_Color color)
    {
        SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);

        SDL_Rect rect = {
            static_cast<int>(box.getCenter().x - box.getSize().x / 2),
            static_cast<int>(box.getCenter().y - box.getSize().y / 2),
            static_cast<int>(box.getSize().x),
            static_cast<int>(box.getSize().y)};
        SDL_RenderDrawRect(renderer, &rect);
    }

    void draw_cross(int x, int y, int size)
    {
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
        // Draw horizontal line
        SDL_RenderDrawLine(renderer, x - size / 2, y, x + size / 2, y);
        // Draw vertical line
        SDL_RenderDrawLine(renderer, x, y - size / 2, x, y + size / 2);
    }

    void draw_text(const std::string &text, int x, int y)
    {
        SDL_Color color = {255, 255, 255, 255};
        SDL_Surface *surface = TTF_RenderText_Solid(font, text.c_str(), color);
        SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);

        int texW = 0;
        int texH = 0;
        SDL_QueryTexture(texture, NULL, NULL, &texW, &texH);
        SDL_Rect dstrect = {x, y, texW, texH};

        SDL_RenderCopy(renderer, texture, NULL, &dstrect);

        SDL_DestroyTexture(texture);
        SDL_FreeSurface(surface);
    }

    void draw_particle(Particle *particle)
    {
        SDL_Color color = {static_cast<Uint8>(particle->color[0] * 255),
                           static_cast<Uint8>(particle->color[1] * 255),
                           static_cast<Uint8>(particle->color[2] * 255),
                           static_cast<Uint8>(particle->color[3] * 255)};
        draw_circle(particle->position.x, particle->position.y, color, particle->radius);
    }

    void init_sdl()
    {
        if (!graphical)
            return;

        std::cout << "Initializing SDL..." << std::endl;
        if (SDL_Init(SDL_INIT_VIDEO) < 0)
        {
            std::cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << std::endl;
            exit(1);
        }

        // ttf init
        if (TTF_Init() < 0)
        {
            std::cerr << "TTF could not initialize! TTF_Error: " << TTF_GetError() << std::endl;
            exit(1);
        }

        window = SDL_CreateWindow("Particle Simulation",
                                  SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                  WINDOW_WIDTH, WINDOW_HEIGHT,
                                  SDL_WINDOW_SHOWN);
        if (!window)
        {
            std::cerr << "Window could not be created! SDL_Error: " << SDL_GetError() << std::endl;
            exit(1);
        }

        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
        if (!renderer)
        {
            std::cerr << "Renderer could not be created! SDL_Error: " << SDL_GetError() << std::endl;
            exit(1);
        }
        SDL_SetHint(SDL_HINT_RENDER_DRIVER, "openGL");

        font = TTF_OpenFont("JetBrainsMono-ExtraLight.ttf", 21);
        if (!font)
        {
            std::cerr << "Error loading font: " << TTF_GetError() << std::endl;
            exit(1);
        }
        std::cout << "SDL initialized" << std::endl;
    }

    void clean_sdl()
    {
        if (!graphical)
            return;

        // cleanup sdl
        TTF_CloseFont(font);
        font = nullptr;
        SDL_DestroyRenderer(renderer);
        renderer = nullptr;
        SDL_DestroyWindow(window);
        window = nullptr;
        TTF_Quit();
        SDL_Quit();
    }
};
